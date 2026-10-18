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
# Low threshold: inlining the linear part creates one traced kernel per
# DISTINCT (group type, L) pair, and problems like LRW1A carry ~70 distinct
# widths of an expensive group type (hours of compilation). Beyond L=4 the
# linear part rides the add_reduce lowering, so all widths of a group type
# collapse into one bucket.
const WIDE_LIN = 4
const ELEM_BUDGET = 40     # max element data fields traced in one kernel

getp(d, s::Symbol) = getproperty(d, s)

# Branch-select multiply (0 * x == 0 even for non-finite x): provided by
# ExaModels since PR #330; a local registration is kept only as a fallback for
# older ExaModels revisions.
@static if isdefined(ExaModels, :bselect)
    const bselect = ExaModels.bselect
else
    bselect(s::Real, x::Real) = ifelse(iszero(s), zero(x), s * x)
    ExaModels.@register_bivariate(bselect, (s, x) -> zero(s), (s, x) -> s,
        (s, x) -> zero(s), (s, x) -> zero(s), (s, x) -> zero(s))
end

# ---------------------------------------------------------------- temp plans

const _NO_FNS = Dict{String,Any}()

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
        # Fortran INTEGER assignment truncates
        nm in et.itemps && haskey(env, nm) && (env[nm] = trunc(env[nm]))
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
function evaltemps!(env, temps, iset::Set{Symbol} = Set{Symbol}())
    for (nm, ent) in temps
        if ent isa CondAssign
            haskey(env, ent.cond) ||
                error("conditional on undefined temporary $(ent.cond)")
            s = tostep(_unint(env[ent.cond]))
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
        # tag integer-typed temporaries so ** lowers to integer power
        nm in iset && !(env[nm] isa IntVal) && (env[nm] = IntVal(env[nm]))
    end
end

"Build the ExaModels node for one element instance at trace time."
function elemnode(x, et::ElType, plan::ElPlan, d, prefix::String, vs::Bool = false,
                  fns::Dict{String,Any} = _NO_FNS)
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
        v = getp(d, Symbol(prefix, "t", k))
        env[nm] = nm in et.itemps ? IntVal(v) : v
    end
    # A caller-supplied Julia function takes precedence over the SIF's own F
    # expression, and is the ONLY way to evaluate a type decoded with
    # external = :defer (its body is opaque Fortran). It receives the env the
    # tracer builds -- elemental and internal variables as ExaModels nodes,
    # element parameters as data fields -- and returns the element's value.
    # Derivatives are ExaModels' problem, so a port supplies the VALUE only,
    # which is why it cannot repeat the hand-coded-derivative errors the
    # Fortran G/H lines are prone to.
    f = get(fns, et.name, nothing)
    f === nothing || return f(env)
    evaltemps!(env, plan.symbolic, et.itemps)
    et.fexpr === nothing && error(
        "element type $(et.name) has no F expression; if the SIF was decoded " *
        "with external = :defer, pass element_fns = Dict(\"$(et.name)\" => f) " *
        "to build_exa")
    return _unint(feval(et.fexpr, env))
end

"Build the ExaModels node for one nontrivial group instance at trace time.
`stages` carries buffered partial sums (linear and/or element chunks)
referenced via the slot field `d.si`; `L` and `ets/plans` describe the parts
traced inline."
function groupnode(x, gt::GType, ets::Vector{ElType}, plans::Vector{ElPlan},
                   L::Int, d, stages = (), vs::Bool = false,
                   fns::Dict{String,Any} = _NO_FNS)
    acc = nothing
    addt(t) = (acc = acc === nothing ? t : acc + t)
    for S in stages
        addt(S[getp(d, :si)])
    end
    for k in 1:L
        addt(getp(d, Symbol("a", k)) * x[getp(d, Symbol("j", k))])
    end
    for (e, et) in enumerate(ets)
        addt(getp(d, Symbol("w", e)) * elemnode(x, et, plans[e], d, "e$(e)_", vs, fns))
    end
    addt(-getp(d, :b))
    env = Dict{Symbol,Any}(gt.gvar => acc)
    for (k, p) in enumerate(gt.gpars)
        env[p] = getp(d, Symbol("g", k))
    end
    evaltemps!(env, gt.temps, gt.itemps)
    return _unint(feval(gt.fexpr, env)) * getp(d, :sinv)
end

# One data row per model term, as a NamedTuple whose fields ExaModels reads
# inside the generator. Float fields are narrowed to the model's element type
# and integer fields left alone: the floats are the SIF's parameters and
# coefficients, which must carry the model's precision, while the integers are
# variable positions, which are positions in any precision. Without the
# conversion a `T = Float32` model holds Float64 parameter tables and computes
# its terms in Float64 before rounding -- on a GPU backend that is Float64
# arithmetic in a Float32 kernel.
function ntvec(::Type{T}, names::Vector{Symbol}, rows::Vector{Vector{Any}}) where {T}
    NT = NamedTuple{Tuple(names)}
    conv(v) = v isa AbstractFloat ? T(v) : v
    return [NT(Tuple(conv(v) for v in r)) for r in rows]
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
    prodflag(prod, c::ExaModels.ExaCore) -> Bool

Resolve the `prod` setting of a model builder for the core `c`.

`prod = true` asks `ExaModels.ExaModel` for the product helper that the
KernelAbstractions `hprod!` / `jprod_nln!` / `jtprod_nln!` overrides dispatch
on. Without that helper those callbacks fall back to the generic methods,
which index the model's arrays elementwise and so fail on a device under
`CUDA.allowscalar(false)`. It is free when `backend === nothing`, where
`ExaModels.build_extension` ignores it.

ExaModels rejects `prod = true` for a core carrying buffered subexpressions,
which wide nontrivial groups emit (see `WIDE_LIN`), so `:auto` requests the
helper exactly where the core admits one. A literal `true` or `false` is
passed through unchanged, so `prod = true` on such a core raises rather than
quietly building a model whose `hprod` does not work on a device.
"""
prodflag(prod, c::ExaModels.ExaCore) = prod === :auto ? isempty(c.subexprs) : prod

"""
    sif_model(path; T = Float64, backend = nothing, prod = :auto, params = Dict(),
              loopfinal = :last) -> ExaModel

Decode a SIF file and build the corresponding ExaModel. Constraint rows appear
in SIF group declaration order (matching SIFDecode/CUTEst ordering).

`prod` is resolved by [`prodflag`](@ref) and forwarded to
`ExaModels.ExaModel`.

Throws [`ExternalFortranError`](@ref) for the 34 collection problems whose
element or group functions are hand-coded Fortran appended to the SIF file:
those are not expressions and cannot be decoded. Evaluating them means
calling the compiled Fortran, which is deliberately NOT this package's
business -- ExaCUTEst.jl handles that tier, so nothing here depends on
CUTEst.jl or on a Fortran toolchain.
"""
function sif_model(path::AbstractString; T = Float64, backend = nothing,
                   prod = :auto, kwargs...)
    return build_exa(sif_decode(path; kwargs...); T, backend, prod)
end

"""
    build_exa(P::SIFProblem; T = Float64, backend = nothing, prod = :auto,
              element_fns = Dict()) -> ExaModel

Build an ExaModel from a decoded problem.

`prod` is resolved by [`prodflag`](@ref) and forwarded to
`ExaModels.ExaModel`.
"""
function build_exa(P::SIFProblem; T = Float64, backend = nothing, prod = :auto,
                   element_fns::AbstractDict = _NO_FNS)
    fns = Dict{String,Any}(String(k) => v for (k, v) in element_fns)
    for (tname, fn) in P.external_types
        startswith(tname, "G#") && error(
            "group type $(tname[3:end]) calls external Fortran ($fn); only " *
            "ELEMENT types can be supplied through element_fns")
        haskey(fns, tname) || error(
            "element type $tname calls external Fortran ($fn) and was decoded " *
            "with external = :defer; supply it as element_fns[\"$tname\"]")
    end
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

    # CUTEst convention: bound magnitudes >= 1e20 are infinite
    for i in 1:ncon
        lcon[i] <= -1e20 && (lcon[i] = T(-Inf))
        ucon[i] >= 1e20 && (ucon[i] = T(Inf))
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

    # --- emit a wide bucket: buffer all elements per type, build each group's
    # argument with one add_reduce (terms are data: O(1) kernel types however
    # many elements/linear terms a group carries), finish with one tiny kernel.
    emitwide = function (sig, bucket, iscon)
        ets, epls, heads, elemrows, lins = bucket
        gt = P.gtypes[sig[2]]
        m = length(heads)
        for (slot, head) in enumerate(heads)
            insert!(head, iscon ? 2 : 1, slot)
        end
        E = length(ets)
        terms = Pair{Int,Any}[]
        if E > 0
            fw = [1 + nelemfields(ets[e], epls[e], hasvs) for e in 1:E]
            offs = cumsum(vcat(0, fw))
            # ONE buffered stage per element TYPE, rows spanning every instance
            # of every group (stages sit in the core's tuple type, so the stage
            # count must not scale with the problem)
            bytype = Dict{String,Vector{Int}}()
            for e in 1:E
                push!(get!(bytype, ets[e].name, Int[]), e)
            end
            for (tname, pos) in bytype
                et = ets[pos[1]]; pl = epls[pos[1]]
                names = Symbol[:w]
                elemfields!(names, et, pl, "", hasvs)
                data = Vector{Any}[]
                owner = Int[]
                for (g, erow) in enumerate(elemrows)
                    for e in pos
                        push!(data, erow[offs[e]+1 : offs[e]+fw[e]])
                        push!(owner, g)
                    end
                end
                local S
                c, S = ExaModels.add_expr(c,
                    Base.Generator(d -> d.w * elemnode(x, et, pl, d, "", hasvs, fns),
                                   ntvec(T, names, data)); buffered = true)
                for (k, g) in enumerate(owner)
                    push!(terms, g => S[k])
                end
            end
        end
        for (g, lin) in enumerate(lins)
            for (a, j) in lin
                push!(terms, g => a * x[j])
            end
        end
        local alpha
        c, alpha = ExaModels.add_reduce(c, m, terms)   # -b applied in the final kernel
        hnames = Symbol[]
        iscon && push!(hnames, :i)
        push!(hnames, :si); push!(hnames, :sinv); push!(hnames, :b)
        for k in 1:length(gt.gpars); push!(hnames, Symbol("g", k)); end
        st = (alpha,)
        noets = ElType[]; nopls = ElPlan[]
        data = ntvec(T, hnames, heads)
        if iscon
            c, _ = ExaModels.add_con!(c, con,
                Base.Generator(d -> d.i => groupnode(x, gt, noets, nopls, 0, d, st, hasvs, fns), data))
        else
            c, _ = ExaModels.add_obj(c,
                Base.Generator(d -> groupnode(x, gt, noets, nopls, 0, d, st, hasvs, fns), data))
        end
    end

    # --- emit constraint terms
    if !isempty(ccst)
        data = ntvec(T, [:i, :k], ccst)
        c, _ = ExaModels.add_con!(c, con, Base.Generator(d -> d.i => d.k, data))
    end
    if !isempty(clin)
        data = ntvec(T, [:i, :j, :a], clin)
        c, _ = ExaModels.add_con!(c, con, Base.Generator(d -> d.i => d.a * x[d.j], data))
    end
    for (tname, rows) in celem
        et = P.eltypes[tname]
        pl = get!(() -> elplan(et), plans, tname)
        names = Symbol[:i, :w]
        elemfields!(names, et, pl, "", hasvs)
        data = ntvec(T, names, rows)
        c, _ = ExaModels.add_con!(c, con,
            Base.Generator(d -> d.i => d.w * elemnode(x, et, pl, d, "", hasvs, fns), data))
    end
    for (sig, (names, ets, epls, L, rows)) in cgrp
        gt = P.gtypes[sig[2]]
        data = ntvec(T, names, rows)
        c, _ = ExaModels.add_con!(c, con,
            Base.Generator(d -> d.i => groupnode(x, gt, ets, epls, L, d, (), hasvs, fns), data))
    end
    for (sig, bucket) in cwide
        emitwide(sig, bucket, true)
    end

    # --- emit objective terms
    if !isempty(olin)
        data = ntvec(T, [:j, :a], olin)
        c, _ = ExaModels.add_obj(c, Base.Generator(d -> d.a * x[d.j], data))
    end
    for (tname, rows) in oelem
        et = P.eltypes[tname]
        pl = get!(() -> elplan(et), plans, tname)
        names = Symbol[:w]
        elemfields!(names, et, pl, "", hasvs)
        data = ntvec(T, names, rows)
        c, _ = ExaModels.add_obj(c,
            Base.Generator(d -> d.w * elemnode(x, et, pl, d, "", hasvs, fns), data))
    end
    for (sig, (names, ets, epls, L, rows)) in ogrp
        gt = P.gtypes[sig[2]]
        data = ntvec(T, names, rows)
        c, _ = ExaModels.add_obj(c,
            Base.Generator(d -> groupnode(x, gt, ets, epls, L, d, (), hasvs, fns), data))
    end
    for (sig, bucket) in owide
        emitwide(sig, bucket, false)
    end
    if !isempty(P.quad)
        qrows = Vector{Any}[Any[i, j,
            VS[i] * VS[j] * ((QUADCONV[] === :half && i == j) ? 0.5 * q : q)] for (i, j, q) in P.quad]
        data = ntvec(T, [:i, :j, :q], qrows)
        c, _ = ExaModels.add_obj(c, Base.Generator(d -> d.q * x[d.i] * x[d.j], data))
    end
    if objconst != 0.0
        c, _ = ExaModels.add_obj(c, Base.Generator(d -> d.k, [(k = T(objconst),)]))
    end

    return ExaModels.ExaModel(c; prod = prodflag(prod, c))
end
