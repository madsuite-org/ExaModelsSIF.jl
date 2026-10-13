# builder.jl — map a SIFProblem IR onto the ExaModels SIMD abstraction.
# Elements/groups are bucketed by type signature; each bucket becomes one
# generator over a concrete Vector{NamedTuple}, so same-type SIF entities
# share one compiled kernel.

# Convention for the QUADRATIC/QUADOBJ section: with :half, an entry (i,j,q)
# contributes q*x[i]*x[j] off-diagonal and 0.5*q*x[i]^2 on the diagonal
# (Hessian-entry reading). Calibrated against SIFDecode via CUTEst.jl.
const QUADCONV = Ref(:half)

# Nontrivial groups wider than WIDE_LIN linear terms route their linear part
# through buffered subexpressions (add_expr(...; buffered = true), ExaModels
# compile-time-prototypes branch): partial sums of CHUNK terms each, evaluated
# once per group into a theta-backed stage, so no traced kernel grows with the
# group's linear width. Measured on this host: L=200 monolithic kernel needs
# hours; chunked+buffered builds in ~16 s with exact derivatives.
const WIDE_LIN = 16
const CHUNK = 16

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

getp(d, s::Symbol) = getproperty(d, s)

"Build the ExaModels node for one element instance at trace time."
function elemnode(x, et::ElType, d, prefix::String)
    env = Dict{Symbol,Any}()
    for (k, ev) in enumerate(et.evars)
        env[ev] = x[getp(d, Symbol(prefix, "v", k))]
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
    for (nm, ast) in et.temps
        env[nm] = feval(ast, env)
    end
    et.fexpr === nothing && error("element type $(et.name) has no F expression")
    return feval(et.fexpr, env)
end

"Build the ExaModels node for one nontrivial group instance at trace time.
`stages` is empty for narrow groups; for wide groups it carries the buffered
partial-sum subexpressions referenced via the slot field `d.si`."
function groupnode(x, gt::GType, ets::Vector{ElType}, L::Int, d, stages = ())
    acc = nothing
    addt(t) = (acc = acc === nothing ? t : acc + t)
    for S in stages
        addt(S[getp(d, :si)])
    end
    for k in 1:L
        addt(getp(d, Symbol("a", k)) * x[getp(d, Symbol("j", k))])
    end
    for (e, et) in enumerate(ets)
        addt(getp(d, Symbol("w", e)) * elemnode(x, et, d, "e$(e)_"))
    end
    addt(-getp(d, :b))
    env = Dict{Symbol,Any}(gt.gvar => acc)
    for (k, p) in enumerate(gt.gpars)
        env[p] = getp(d, Symbol("g", k))
    end
    for (nm, ast) in gt.temps
        env[nm] = feval(ast, env)
    end
    return feval(gt.fexpr, env) * getp(d, :sinv)
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

"Element data row fields appended to `row` for element `ename`; returns ElType."
function pushelem!(row, P::SIFProblem, ename::String)
    e = get(P.elements, ename, nothing)
    e === nothing && error("group references undefined element $ename")
    et = eltypeof(P, e)
    for ev in et.evars
        haskey(e.varmap, ev) || error("element $(e.name): variable $ev unassigned")
        push!(row, e.varmap[ev])
    end
    for p in et.epars
        haskey(e.parvals, p) || error("element $(e.name): parameter $p unassigned")
        push!(row, e.parvals[p])
    end
    return et
end

function elemfields!(names, et::ElType, prefix::String)
    for k in 1:length(et.evars)
        push!(names, Symbol(prefix, "v", k))
    end
    for k in 1:length(et.epars)
        push!(names, Symbol(prefix, "p", k))
    end
end

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

    connames = [gn for gn in P.grouporder if P.groups[gn].kind != 'N']
    rowof = Dict(gn => i for (i, gn) in enumerate(connames))
    ncon = length(connames)

    # --- constraint bounds
    lcon = zeros(T, ncon); ucon = zeros(T, ncon)
    for (i, gn) in enumerate(connames)
        g = P.groups[gn]
        triv = gtypename(P, g) === nothing
        s = g.scale
        triv && s <= 0.0 && error("non-positive scale on trivial constraint group $gn")
        # SIFDecode/CUTEst convention: the group constant is folded into the
        # constraint value (c(x) = alpha - b), so bounds are built around 0.
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
    celem = Dict{String,Vector{Vector{Any}}}()   # etype => rows (i, w, v.., p..)
    oelem = Dict{String,Vector{Vector{Any}}}()   # etype => rows (w, v.., p..)
    cgrp = Dict{Any,Any}()                    # narrow signature => (names, ets, L, rows)
    ogrp = Dict{Any,Any}()
    cwide = Dict{Any,Any}()                   # wide signature => (names, ets, rows, lins)
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
                iscon ? push!(clin, Any[rowof[gn], j, a / s]) :
                        push!(olin, Any[j, a / s])
            end
            for (ename, w) in g.elems
                e = P.elements[ename]
                et = eltypeof(P, e)
                row = Any[]
                iscon && push!(row, rowof[gn])
                push!(row, w / s)
                pushelem!(row, P, ename)
                push!(get!(iscon ? celem : oelem, et.name, Vector{Any}[]), row)
            end
            if iscon
                g.constant != 0.0 && push!(ccst, Any[rowof[gn], -g.constant / s])
            else
                objconst -= g.constant / s
            end
        else
            gt = P.gtypes[tn]
            etnames = Tuple(eltypeof(P, P.elements[en]).name for (en, _) in g.elems)
            nelemfields = sum(Int[length(eltypeof(P, P.elements[en]).evars) +
                length(eltypeof(P, P.elements[en]).epars) for (en, _) in g.elems])
            wide = length(lin) > WIDE_LIN
            # Even with the linear part buffered, the element part is traced in
            # one kernel; refuse cleanly when that alone is too wide.
            nelemfields + (wide ? cld(length(lin), CHUNK) : length(lin)) > 64 &&
                error("unsupported: nontrivial group $gn too wide for " *
                      "single-kernel tracing (elem fields=$nelemfields, nlin=$(length(lin)))")
            row = Any[]
            iscon && push!(row, rowof[gn])
            push!(row, 1.0 / s)
            push!(row, g.constant)
            for p in gt.gpars
                haskey(g.gparvals, p) || error("group $gn: group parameter $p unassigned")
                push!(row, g.gparvals[p])
            end
            if wide
                ets = ElType[]
                for (ename, w) in g.elems
                    push!(row, w)
                    push!(ets, pushelem!(row, P, ename))
                end
                sig = (iscon, tn, etnames)
                if !haskey(iscon ? cwide : owide, sig)
                    names = Symbol[]
                    iscon && push!(names, :i)
                    push!(names, :si); push!(names, :sinv); push!(names, :b)
                    for k in 1:length(gt.gpars); push!(names, Symbol("g", k)); end
                    for (e, et) in enumerate(ets)
                        push!(names, Symbol("w", e))
                        elemfields!(names, et, "e$(e)_")
                    end
                    (iscon ? cwide : owide)[sig] = (names, ets, Vector{Any}[], Vector{Vector{Tuple{Float64,Int}}}())
                end
                b = (iscon ? cwide : owide)[sig]
                push!(row, 0)  # placeholder; si assigned at emit (slot order)
                # move si placeholder into position: rebuild row with si after i
                push!(b[3], row)
                push!(b[4], [(a, j) for (j, a) in lin])
            else
                for (k, (j, a)) in enumerate(lin)
                    push!(row, a); push!(row, j)
                end
                ets = ElType[]
                for (e, (ename, w)) in enumerate(g.elems)
                    push!(row, w)
                    push!(ets, pushelem!(row, P, ename))
                end
                sig = (iscon, tn, length(lin), etnames)
                if !haskey(iscon ? cgrp : ogrp, sig)
                    names = Symbol[]
                    iscon && push!(names, :i)
                    push!(names, :sinv); push!(names, :b)
                    for k in 1:length(gt.gpars); push!(names, Symbol("g", k)); end
                    for k in 1:length(lin); push!(names, Symbol("a", k)); push!(names, Symbol("j", k)); end
                    for (e, et) in enumerate(ets)
                        push!(names, Symbol("w", e))
                        elemfields!(names, et, "e$(e)_")
                    end
                    (iscon ? cgrp : ogrp)[sig] = (names, ets, length(lin), Vector{Any}[])
                end
                push!((iscon ? cgrp : ogrp)[sig][4], row)
            end
        end
    end

    # --- emit a wide bucket: buffered partial-sum stages + one narrow kernel
    emitwide = function (sig, bucket, iscon)
        names, ets, rows, lins = bucket
        gt = P.gtypes[sig[2]]
        m = length(rows)
        # assign slots in bucket order; si sits right after the optional i
        for (slot, row) in enumerate(rows)
            row[end] = 0  # drop placeholder tail
            pop!(row)
            insert!(row, iscon ? 2 : 1, slot)
        end
        Lmax = maximum(length.(lins))
        nlay = cld(Lmax, CHUNK)
        ker = LinChunk(x, CHUNK)
        cnames = vcat([Symbol("a", k) for k in 1:CHUNK], [Symbol("j", k) for k in 1:CHUNK])
        CNT = NamedTuple{Tuple(cnames)}
        stages = Any[]
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
        st = Tuple(stages)
        data = ntvec(names, rows)
        if iscon
            c, _ = ExaModels.add_con!(c, con,
                Base.Generator(d -> d.i => groupnode(x, gt, ets, 0, d, st), data))
        else
            c, _ = ExaModels.add_obj(c,
                Base.Generator(d -> groupnode(x, gt, ets, 0, d, st), data))
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
        names = Symbol[:i, :w]
        elemfields!(names, et, "")
        data = ntvec(names, rows)
        c, _ = ExaModels.add_con!(c, con,
            Base.Generator(d -> d.i => d.w * elemnode(x, et, d, ""), data))
    end
    for (sig, (names, ets, L, rows)) in cgrp
        gt = P.gtypes[sig[2]]
        data = ntvec(names, rows)
        c, _ = ExaModels.add_con!(c, con,
            Base.Generator(d -> d.i => groupnode(x, gt, ets, L, d), data))
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
        names = Symbol[:w]
        elemfields!(names, et, "")
        data = ntvec(names, rows)
        c, _ = ExaModels.add_obj(c,
            Base.Generator(d -> d.w * elemnode(x, et, d, ""), data))
    end
    for (sig, (names, ets, L, rows)) in ogrp
        gt = P.gtypes[sig[2]]
        data = ntvec(names, rows)
        c, _ = ExaModels.add_obj(c,
            Base.Generator(d -> groupnode(x, gt, ets, L, d), data))
    end
    for (sig, bucket) in owide
        emitwide(sig, bucket, false)
    end
    if !isempty(P.quad)
        qrows = Vector{Any}[Any[i, j,
            (QUADCONV[] === :half && i == j) ? 0.5 * q : q] for (i, j, q) in P.quad]
        data = ntvec([:i, :j, :q], qrows)
        c, _ = ExaModels.add_obj(c, Base.Generator(d -> d.q * x[d.i] * x[d.j], data))
    end
    if objconst != 0.0
        c, _ = ExaModels.add_obj(c, Base.Generator(d -> d.k, [(k = T(objconst),)]))
    end

    return ExaModels.ExaModel(c)
end
