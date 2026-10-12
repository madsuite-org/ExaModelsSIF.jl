# builder.jl — map a SIFProblem IR onto the ExaModels SIMD abstraction.
# Elements/groups are bucketed by type signature; each bucket becomes one
# generator over a concrete Vector{NamedTuple}, so same-type SIF entities
# share one compiled kernel.

# Convention for the QUADRATIC/QUADOBJ section: with :half, an entry (i,j,q)
# contributes q*x[i]*x[j] off-diagonal and 0.5*q*x[i]^2 on the diagonal
# (Hessian-entry reading). Calibrated against SIFDecode via CUTEst.jl.
const QUADCONV = Ref(:half)

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

"Build the ExaModels node for one nontrivial group instance at trace time."
function groupnode(x, gt::GType, ets::Vector{ElType}, L::Int, d)
    acc = nothing
    addt(t) = (acc = acc === nothing ? t : acc + t)
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
    cgrp = Dict{Any,Any}()                    # signature => (names, ets, L, rows)
    ogrp = Dict{Any,Any}()
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
            sig = (iscon, tn, length(lin), etnames)
            row = Any[]
            iscon && push!(row, rowof[gn])
            push!(row, 1.0 / s)
            push!(row, g.constant)
            for p in gt.gpars
                haskey(g.gparvals, p) || error("group $gn: group parameter $p unassigned")
                push!(row, g.gparvals[p])
            end
            for (k, (j, a)) in enumerate(lin)
                push!(row, a); push!(row, j)
            end
            ets = ElType[]
            for (e, (ename, w)) in enumerate(g.elems)
                push!(row, w)
                push!(ets, pushelem!(row, P, ename))
            end
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
