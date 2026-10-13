# builder.jl — map a SIFProblem IR onto the ExaModels SIMD abstraction.
# Elements/groups are bucketed by type signature; each bucket becomes one
# generator over a concrete Vector{NamedTuple}, so same-type SIF entities
# share one compiled kernel.

# Convention for the QUADRATIC/QUADOBJ section: with :half, an entry (i,j,q)
# contributes q*x[i]*x[j] off-diagonal and 0.5*q*x[i]^2 on the diagonal
# (Hessian-entry reading). Calibrated against SIFDecode via CUTEst.jl.
const QUADCONV = Ref(:half)

# Nontrivial groups wider than WIDE_LIN linear terms (or ELEM_BUDGET element
# fields) route the wide part through buffered subexpressions
# (add_expr(...; buffered = true), ExaModels compile-time-prototypes branch):
# partial sums evaluated once per group into a theta-backed stage, so no traced
# kernel grows with the group's width. Measured on 002: L=200 monolithic kernel
# needs hours; chunked+buffered builds in ~16 s with exact derivatives.
const WIDE_LIN = 16
const CHUNK = 16
const ELEM_BUDGET = 40     # max element data fields traced in one kernel
const ELEM_CHUNK = 24      # target element-field count per buffered elem stage

getp(d, s::Symbol) = getproperty(d, s)

# Branch-select multiply: bselect(s, x) = s*x with the convention 0*x = 0 even
# for x = Inf/NaN, so an INACTIVE conditional branch cannot poison the blend
# (QC-class guarded logs). All derivative stubs are inert: the step factor s
# carries zero derivative by construction, and d/dx is s itself.
bselect(s, x) = iszero(s) ? zero(x) : s * x
_bs_d1(s, x) = zero(s)
_bs_d2(s, x) = s
_bs_d11(s, x) = zero(s)
_bs_d12(s, x) = zero(s)
_bs_d22(s, x) = zero(s)
ExaModels.@register_bivariate(bselect, _bs_d1, _bs_d2, _bs_d11, _bs_d12, _bs_d22)

struct LinChunk{X}
    x::X
    C::Int
end
function (k::LinChunk)(d)
    acc = getproperty(d, :a1) * k.x[getproperty(d, :j1)]
    for q in 2:k.C
        acc = acc + getproperty(d, Symbol("a", q)) * k.x[getproperty(d, Symbol("j", q))]
    end
    return acc
end

# ---------------------------------------------------------------- temp plans

"Evaluation plan for an element type's temporaries: entries computable from
element parameters alone are evaluated numerically per instance at build time
(enabling parameter-only conditionals); the rest trace symbolically."
struct ElPlan
    paramonly::Vector{Pair{Symbol,Any}}
    inject::Vector{Symbol}     # param-only temps consumed by F / symbolic temps
    symbolic::Vector{Pair{Symbol,Any}}
end

function elplan(et::ElType)
    pset = Set(et.epars)
    ponames = Set{Symbol}()
    paramonly = Pair{Symbol,Any}[]
    symbolic = Pair{Symbol,Any}[]
    for (nm, ent) in et.temps
        if all(d -> d in pset || d in ponames, entdeps(ent))
            push!(paramonly, nm => ent)
            push!(ponames, nm)
        else
            push!(symbolic, nm => ent)
        end
    end
    used = et.fexpr === nothing ? Set{Symbol}() : freesyms(et.fexpr)
    for (nm, ent) in symbolic
        union!(used, entdeps(ent))
    end
    inject = Symbol[]
    for (nm, _) in paramonly
        nm in used && !(nm in inject) && push!(inject, nm)
    end
    return ElPlan(paramonly, inject, symbolic)
end

"Numerically evaluate the param-only temporaries of one element instance;
returns the injected field values in plan.inject order."
function evalparamtemps(plan::ElPlan, et::ElType, e::SIFElement)
    env = Dict{Symbol,Any}()
    for p in et.epars
        haskey(e.parvals, p) || error("element $(e.name): parameter $p unassigned")
        env[p] = e.parvals[p]
    end
    for (nm, ent) in plan.paramonly
        if ent isa CondAssign
            haskey(env, ent.cond) ||
                error("element $(e.name): conditional on undefined temp $(ent.cond)")
            (env[ent.cond] == ent.sense) && (env[nm] = feval(ent.ast, env))
        else
            env[nm] = feval(ent, env)
        end
    end
    vals = Float64[]
    for nm in plan.inject
        haskey(env, nm) ||
            error("element $(e.name): param-only temp $nm never assigned")
        push!(vals, Float64(env[nm]))
    end
    return vals
end

"Evaluate temporaries symbolically; conditional assignments blend branches
with exact step factors (see relstep), so variable-dependent conditionals
reproduce SIFDecode's branch-selected values and derivatives."
function evaltemps!(env, temps)
    for (nm, ent) in temps
        if ent isa CondAssign
            haskey(env, ent.cond) ||
                error("conditional on undefined temporary $(ent.cond)")
            s = tostep(env[ent.cond])
            v = feval(ent.ast, env)
            w = ent.sense ? s : (s isa Number ? 1.0 - s : 1.0 - s)
            if s isa Number
                (w == 1.0) && (env[nm] = v)
                (w == 0.0) || (w == 1.0) ||
                    error("non-binary numeric condition for $(ent.cond)")
            else
                prev = get(env, nm, nothing)
                if prev === nothing
                    env[nm] = bselect(w, v)
                else
                    env[nm] = bselect(w, v) + bselect(1.0 - w, prev)
                end
            end
        else
            env[nm] = feval(ent, env)
        end
    end
end

"Build the ExaModels node for one element instance at trace time."
function elemnode(x, et::ElType, plan::ElPlan, d, prefix::String, vs::Bool = false)
    env = Dict{Symbol,Any}()
    for (k, ev) in enumerate(et.evars)
        xn = x[getp(d, Symbol(prefix, "v", k))]
        env[ev] = vs ? getp(d, Symbol(prefix, "s", k)) * xn : xn
    end
    for iv in et.ivars
        rows = get(et.rrows, iv, nothing)
        rows === nothing && error("internal variable $iv of $(et.name) lacks R rows")
        acc = nothing
        for (ev, cf) in rows
            t = cf == 1.0 ? env[ev] : cf * env[ev]
            acc = acc === nothing ? t : acc + t
        end
        env[iv] = acc
    end
    for (k, p) in enumerate(et.epars)
        env[p] = getp(d, Symbol(prefix, "p", k))
    end
    for (k, nm) in enumerate(plan.inject)
        env[nm] = getp(d, Symbol(prefix, "t", k))
    end
    evaltemps!(env, plan.symbolic)
    et.fexpr === nothing && error("element type $(et.name) has no F expression")
    return feval(et.fexpr, env)
end

"Build the ExaModels node for one nontrivial group instance at trace time.
`stages` carries buffered partial sums (linear and/or element chunks)
referenced via the slot field `d.si`; `L` and `ets/plans` describe the parts
traced inline."
function groupnode(x, gt::GType, ets::Vector{ElType}, plans::Vector{ElPlan},
                   L::Int, d, stages = (), vs::Bool = false)
    acc = nothing
    addt(t) = (acc = acc === nothing ? t : acc + t)
    for S in stages
        addt(S[getp(d, :si)])
    end
    for k in 1:L
        addt(getp(d, Symbol("a", k)) * x[getp(d, Symbol("j", k))])
    end
    for (e, et) in enumerate(ets)
        addt(getp(d, Symbol("w", e)) * elemnode(x, et, plans[e], d, "e$(e)_", vs))
    end
    addt(-getp(d, :b))
    env = Dict{Symbol,Any}(gt.gvar => acc)
    for (k, p) in enumerate(gt.gpars)
        env[p] = getp(d, Symbol("g", k))
    end
    evaltemps!(env, gt.temps)
    return feval(gt.fexpr, env) * getp(d, :sinv)
end

"Kernel for one buffered element-chunk stage: sum of the chunk's weighted
elements. Field names are chunk-RELATIVE (c1_, c2_, ...), so chunks with the
same element-type composition share one traced tree type and compile once."
struct ElemChunk{X}
    x::X
    ets::Vector{ElType}
    plans::Vector{ElPlan}
    vs::Bool
end
function (k::ElemChunk)(d)
    acc = nothing
    for rel in 1:length(k.ets)
        t = getp(d, Symbol("w", rel)) * elemnode(k.x, k.ets[rel], k.plans[rel], d, "c$(rel)_", k.vs)
        acc = acc === nothing ? t : acc + t
    end
    return acc
end

function ntvec(names::Vector{Symbol}, rows::Vector{Vector{Any}})
    NT = NamedTuple{Tuple(names)}
    return [NT(Tuple(r)) for r in rows]
end

function eltypeof(P::SIFProblem, e::SIFElement)
    tn = isempty(e.etype) ? P.default_etype : e.etype
    tn === nothing && error("element $(e.name) has no type and no default exists")
    haskey(P.eltypes, tn) || error("element type $tn has no function body")
    return P.eltypes[tn]
end

function gtypename(P::SIFProblem, g::SIFGroup)
    t = g.gtype !== nothing ? g.gtype : P.default_gtype
    (t === nothing || t == "TRIVIAL") && return nothing
    haskey(P.gtypes, t) || error("group type $t has no function body")
    return t
end

"Append one element instance's data fields (variable indices, parameters,
injected temp values) to `row`; returns (ElType, ElPlan, nfields)."
function pushelem!(row, P::SIFProblem, plans::Dict{String,ElPlan}, ename::String, vs::Bool)
    e = get(P.elements, ename, nothing)
    e === nothing && error("group references undefined element $ename")
    et = eltypeof(P, e)
    plan = get!(() -> elplan(et), plans, et.name)
    for ev in et.evars
        haskey(e.varmap, ev) || error("element $(e.name): variable $ev unassigned")
        push!(row, e.varmap[ev])
    end
    for p in et.epars
        haskey(e.parvals, p) || error("element $(e.name): parameter $p unassigned")
        push!(row, e.parvals[p])
    end
    tv = evalparamtemps(plan, et, e)
    append!(row, tv)
    if vs
        for ev in et.evars
            push!(row, P.varscale[e.varmap[ev]])
        end
    end
    return et, plan, nelemfields(et, plan, vs)
end

function elemfields!(names, et::ElType, plan::ElPlan, prefix::String, vs::Bool)
    for k in 1:length(et.evars)
        push!(names, Symbol(prefix, "v", k))
    end
    for k in 1:length(et.epars)
        push!(names, Symbol(prefix, "p", k))
    end
    for k in 1:length(plan.inject)
        push!(names, Symbol(prefix, "t", k))
    end
    if vs
        for k in 1:length(et.evars)
            push!(names, Symbol(prefix, "s", k))
        end
    end
end

nelemfields(et::ElType, plan::ElPlan, vs::Bool) =
    (vs ? 2 : 1) * length(et.evars) + length(et.epars) + length(plan.inject)

combinelin(g::SIFGroup) = begin
    acc = Dict{Int,Float64}()
    order = Int[]
    for (j, a) in g.lin
        haskey(acc, j) || push!(order, j)
        acc[j] = get(acc, j, 0.0) + a
    end
    [(j, acc[j]) for j in order]
end

"""
    sif_model(path; T = Float64, backend = nothing, loopfinal = :last) -> ExaModel

Decode a SIF file and build the corresponding ExaModel. Constraint rows appear
in SIF group declaration order (matching SIFDecode/CUTEst ordering).
"""
function sif_model(path::AbstractString; T = Float64, backend = nothing, kwargs...)
    P = sif_decode(path; kwargs...)
    return build_exa(P; T, backend)
end

function build_exa(P::SIFProblem; T = Float64, backend = nothing)
    nvar = length(P.varnames)
    c = ExaModels.ExaCore(T; backend)
    c, x = ExaModels.add_var(c, nvar; start = P.x0, lvar = P.lvar, uvar = P.uvar)
    plans = Dict{String,ElPlan}()
    VS = P.varscale
    hasvs = any(!=(1.0), VS)

    connames = [gn for gn in P.grouporder if P.groups[gn].kind != 'N']
    rowof = Dict(gn => i for (i, gn) in enumerate(connames))
    ncon = length(connames)

    # --- constraint bounds (constant folded into the constraint value)
    lcon = zeros(T, ncon); ucon = zeros(T, ncon)
    for (i, gn) in enumerate(connames)
        g = P.groups[gn]
        triv = gtypename(P, g) === nothing
        s = g.scale
        triv && s <= 0.0 && error("non-positive scale on trivial constraint group $gn")
        bc = 0.0
        r = g.range_ === nothing ? nothing : (triv ? g.range_ / s : g.range_)
        if g.kind == 'E'
            if r === nothing
                lcon[i] = bc; ucon[i] = bc
            elseif r >= 0
                lcon[i] = bc; ucon[i] = bc + r
            else
                lcon[i] = bc + r; ucon[i] = bc
            end
        elseif g.kind == 'L'
            lcon[i] = r === nothing ? T(-Inf) : bc - abs(r); ucon[i] = bc
        elseif g.kind == 'G'
            lcon[i] = bc; ucon[i] = r === nothing ? T(Inf) : bc + abs(r)
        else
            error("unexpected constraint kind $(g.kind)")
        end
    end

    con = nothing
    if ncon > 0
        c, con = ExaModels.add_con(c, ncon; lcon = lcon, ucon = ucon)
    end

    # --- buckets
    clin = Vector{Any}[]                      # (i, j, a)
    ccst = Vector{Any}[]                      # (i, k) constant term per row
    olin = Vector{Any}[]                      # (j, a)
    celem = Dict{String,Vector{Vector{Any}}}()   # etype => rows (i, w, fields..)
    oelem = Dict{String,Vector{Vector{Any}}}()   # etype => rows (w, fields..)
    cgrp = Dict{Any,Any}()                    # narrow sig => (names, ets, plans, L, rows)
    ogrp = Dict{Any,Any}()
    cwide = Dict{Any,Any}()                   # wide sig => (ets, plans, heads, elemrows, lins)
    owide = Dict{Any,Any}()
    objconst = 0.0

    for gn in P.grouporder
        g = P.groups[gn]
        iscon = g.kind != 'N'
        s = g.scale
        lin = combinelin(g)
        tn = gtypename(P, g)
        if tn === nothing
            for (j, a) in lin
                iscon ? push!(clin, Any[rowof[gn], j, a * VS[j] / s]) :
                        push!(olin, Any[j, a * VS[j] / s])
            end
            for (ename, w) in g.elems
                e = P.elements[ename]
                et0 = eltypeof(P, e)
                row = Any[]
                iscon && push!(row, rowof[gn])
                push!(row, w / s)
                pushelem!(row, P, plans, ename, hasvs)
                push!(get!(iscon ? celem : oelem, et0.name, Vector{Any}[]), row)
            end
            if iscon
                g.constant != 0.0 && push!(ccst, Any[rowof[gn], -g.constant / s])
            else
                objconst -= g.constant / s
            end
        else
            gt = P.gtypes[tn]
            etnames = Tuple(eltypeof(P, P.elements[en]).name for (en, _) in g.elems)
            # field counts decide narrow vs wide lowering
            efields = Int[]
            for (en, _) in g.elems
                e = P.elements[en]
                et0 = eltypeof(P, e)
                pl = get!(() -> elplan(et0), plans, et0.name)
                push!(efields, 1 + nelemfields(et0, pl, hasvs))
            end
            wide = length(lin) > WIDE_LIN || sum(efields; init = 0) > ELEM_BUDGET
            if wide
                sig = (iscon, tn, etnames)
                head = Any[]
                iscon && push!(head, rowof[gn])
                push!(head, 1.0 / s)
                push!(head, g.constant)
                for p in gt.gpars
                    haskey(g.gparvals, p) || error("group $gn: group parameter $p unassigned")
                    push!(head, g.gparvals[p])
                end
                erow = Any[]
                ets = ElType[]; epls = ElPlan[]
                for (ename, w) in g.elems
                    push!(erow, w)
                    et0, pl, _ = pushelem!(erow, P, plans, ename, hasvs)
                    push!(ets, et0); push!(epls, pl)
                end
                b = get!(() -> (ets, epls, Vector{Any}[], Vector{Any}[],
                                Vector{Vector{Tuple{Float64,Int}}}()), iscon ? cwide : owide, sig)
                push!(b[3], head)
                push!(b[4], erow)
                push!(b[5], [(a * VS[j], j) for (j, a) in lin])
            else
                row = Any[]
                iscon && push!(row, rowof[gn])
                push!(row, 1.0 / s)
                push!(row, g.constant)
                for p in gt.gpars
                    haskey(g.gparvals, p) || error("group $gn: group parameter $p unassigned")
                    push!(row, g.gparvals[p])
                end
                for (j, a) in lin
                    push!(row, a * VS[j]); push!(row, j)
                end
                ets = ElType[]; epls = ElPlan[]
                for (ename, w) in g.elems
                    push!(row, w)
                    et0, pl, _ = pushelem!(row, P, plans, ename, hasvs)
                    push!(ets, et0); push!(epls, pl)
                end
                sig = (iscon, tn, length(lin), etnames)
                if !haskey(iscon ? cgrp : ogrp, sig)
                    names = Symbol[]
                    iscon && push!(names, :i)
                    push!(names, :sinv); push!(names, :b)
                    for k in 1:length(gt.gpars); push!(names, Symbol("g", k)); end
                    for k in 1:length(lin); push!(names, Symbol("a", k)); push!(names, Symbol("j", k)); end
                    for (e, et0) in enumerate(ets)
                        push!(names, Symbol("w", e))
                        elemfields!(names, et0, epls[e], "e$(e)_", hasvs)
                    end
                    (iscon ? cgrp : ogrp)[sig] = (names, ets, epls, length(lin), Vector{Any}[])
                end
                push!((iscon ? cgrp : ogrp)[sig][5], row)
            end
        end
    end

    # --- emit a wide bucket: buffered partial-sum stages + one narrow kernel
    emitwide = function (sig, bucket, iscon)
        ets, epls, heads, elemrows, lins = bucket
        gt = P.gtypes[sig[2]]
        m = length(heads)
        for (slot, head) in enumerate(heads)
            insert!(head, iscon ? 2 : 1, slot)
        end
        stages = Any[]
        # linear stages
        Lmax = isempty(lins) ? 0 : maximum(length.(lins))
        if Lmax > 0
            nlay = cld(Lmax, CHUNK)
            ker = LinChunk(x, CHUNK)
            cnames = vcat([Symbol("a", k) for k in 1:CHUNK], [Symbol("j", k) for k in 1:CHUNK])
            CNT = NamedTuple{Tuple(cnames)}
            for t in 1:nlay
                layer = Vector{Any}[]
                for terms in lins
                    ch = terms[min((t-1)*CHUNK+1, length(terms)+1) : min(t*CHUNK, length(terms))]
                    chv = Any[]
                    for (a, j) in ch; push!(chv, a); end
                    for _ in 1:(CHUNK - length(ch)); push!(chv, 0.0); end
                    for (a, j) in ch; push!(chv, j); end
                    for _ in 1:(CHUNK - length(ch)); push!(chv, 1); end
                    push!(layer, chv)
                end
                data = [CNT(Tuple(r)) for r in layer]
                local S
                c, S = ExaModels.add_expr(c, Base.Generator(ker, data); buffered = true)
                push!(stages, S)
            end
        end
        # element stages: partition element positions into field-budget chunks
        E = length(ets)
        if E > 0
            fw = [1 + nelemfields(ets[e], epls[e], hasvs) for e in 1:E]
            chunks = Vector{Int}[]
            cur = Int[]; acc = 0
            for e in 1:E
                if !isempty(cur) && acc + fw[e] > ELEM_CHUNK
                    push!(chunks, cur); cur = Int[]; acc = 0
                end
                push!(cur, e); acc += fw[e]
            end
            isempty(cur) || push!(chunks, cur)
            # offsets of each element's fields within an elemrow
            offs = cumsum(vcat(0, fw))
            for pos in chunks
                names = Symbol[]
                for (rel, e) in enumerate(pos)
                    push!(names, Symbol("w", rel))
                    elemfields!(names, ets[e], epls[e], "c$(rel)_", hasvs)
                end
                data = Vector{Any}[]
                for erow in elemrows
                    rowv = Any[]
                    for e in pos
                        append!(rowv, erow[offs[e]+1 : offs[e]+fw[e]])
                    end
                    push!(data, rowv)
                end
                ker = ElemChunk(x, ets[pos], epls[pos], hasvs)
                local S
                c, S = ExaModels.add_expr(c, Base.Generator(ker, ntvec(names, data)); buffered = true)
                push!(stages, S)
            end
        end
        length(stages) > 64 && error("unsupported: wide group needs " *
            "$(length(stages)) buffered stages (reduction-tree lowering pending)")
        hnames = Symbol[]
        iscon && push!(hnames, :i)
        push!(hnames, :si); push!(hnames, :sinv); push!(hnames, :b)
        for k in 1:length(gt.gpars); push!(hnames, Symbol("g", k)); end
        st = Tuple(stages)
        noets = ElType[]; nopls = ElPlan[]
        data = ntvec(hnames, heads)
        if iscon
            c, _ = ExaModels.add_con!(c, con,
                Base.Generator(d -> d.i => groupnode(x, gt, noets, nopls, 0, d, st, hasvs), data))
        else
            c, _ = ExaModels.add_obj(c,
                Base.Generator(d -> groupnode(x, gt, noets, nopls, 0, d, st, hasvs), data))
        end
    end

    # --- emit constraint terms
    if !isempty(ccst)
        data = ntvec([:i, :k], ccst)
        c, _ = ExaModels.add_con!(c, con, Base.Generator(d -> d.i => d.k, data))
    end
    if !isempty(clin)
        data = ntvec([:i, :j, :a], clin)
        c, _ = ExaModels.add_con!(c, con, Base.Generator(d -> d.i => d.a * x[d.j], data))
    end
    for (tname, rows) in celem
        et = P.eltypes[tname]
        pl = get!(() -> elplan(et), plans, tname)
        names = Symbol[:i, :w]
        elemfields!(names, et, pl, "", hasvs)
        data = ntvec(names, rows)
        c, _ = ExaModels.add_con!(c, con,
            Base.Generator(d -> d.i => d.w * elemnode(x, et, pl, d, "", hasvs), data))
    end
    for (sig, (names, ets, epls, L, rows)) in cgrp
        gt = P.gtypes[sig[2]]
        data = ntvec(names, rows)
        c, _ = ExaModels.add_con!(c, con,
            Base.Generator(d -> d.i => groupnode(x, gt, ets, epls, L, d, (), hasvs), data))
    end
    for (sig, bucket) in cwide
        emitwide(sig, bucket, true)
    end

    # --- emit objective terms
    if !isempty(olin)
        data = ntvec([:j, :a], olin)
        c, _ = ExaModels.add_obj(c, Base.Generator(d -> d.a * x[d.j], data))
    end
    for (tname, rows) in oelem
        et = P.eltypes[tname]
        pl = get!(() -> elplan(et), plans, tname)
        names = Symbol[:w]
        elemfields!(names, et, pl, "", hasvs)
        data = ntvec(names, rows)
        c, _ = ExaModels.add_obj(c,
            Base.Generator(d -> d.w * elemnode(x, et, pl, d, "", hasvs), data))
    end
    for (sig, (names, ets, epls, L, rows)) in ogrp
        gt = P.gtypes[sig[2]]
        data = ntvec(names, rows)
        c, _ = ExaModels.add_obj(c,
            Base.Generator(d -> groupnode(x, gt, ets, epls, L, d, (), hasvs), data))
    end
    for (sig, bucket) in owide
        emitwide(sig, bucket, false)
    end
    if !isempty(P.quad)
        qrows = Vector{Any}[Any[i, j,
            VS[i] * VS[j] * ((QUADCONV[] === :half && i == j) ? 0.5 * q : q)] for (i, j, q) in P.quad]
        data = ntvec([:i, :j, :q], qrows)
        c, _ = ExaModels.add_obj(c, Base.Generator(d -> d.q * x[d.i] * x[d.j], data))
    end
    if objconst != 0.0
        c, _ = ExaModels.add_obj(c, Base.Generator(d -> d.k, [(k = T(objconst),)]))
    end

    return ExaModels.ExaModel(c)
end
