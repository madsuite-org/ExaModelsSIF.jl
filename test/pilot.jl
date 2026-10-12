# test/pilot.jl — M1 pilot validation: compare every decoded model against the
# Fortran reference (CUTEst.jl) on metadata, values, and derivatives.
using ExaModelsSIF, ExaModels, CUTEst, NLPModels
using Random, Printf, LinearAlgebra, SparseArrays

ENV["MASTSIF"] = get(ENV, "MASTSIF", expanduser("~/git/pkg/cutest-sif"))

const PILOT = ["ROSENBR", "HS1", "HS2", "HS3", "HS4", "HS5", "HS6", "HS7",
    "HS10", "HS11", "HS14", "HS21", "HS26", "HS28", "HS35", "HS45",
    "BOOTH", "ZANGWIL3", "TAME", "ALSOTAME", "DEGDIAG", "GENHS28"]

relerr(a, b) = norm(a - b, Inf) / (1.0 + max(norm(a, Inf), norm(b, Inf)))

function spjac(nlp, x)
    r, c = jac_structure(nlp)
    v = jac_coord(nlp, x)
    sparse(r, c, v, nlp.meta.ncon, nlp.meta.nvar)
end

function sphess(nlp, x, y)
    r, c = hess_structure(nlp)
    v = hess_coord(nlp, x, y)
    # canonical lower triangle
    rl = max.(r, c); cl = min.(r, c)
    sparse(rl, cl, v, nlp.meta.nvar, nlp.meta.nvar)
end

# CUTEst default ordering: linear constraints first, equations first within
# each class (efirst = lfirst = true). Returns perm with natural-order index
# of each reference row.
function cutest_perm(P::ExaModelsSIF.SIFProblem)
    connames = [gn for gn in P.grouporder if P.groups[gn].kind != 'N']
    islin = Bool[]; iseq = Bool[]
    for gn in connames
        g = P.groups[gn]
        tn = ExaModelsSIF.gtypename(P, g)
        push!(islin, tn === nothing && isempty(g.elems))
        push!(iseq, g.kind == 'E')
    end
    n = length(connames)
    order = Int[]
    for (l, e) in ((true, true), (true, false), (false, true), (false, false))
        append!(order, [i for i in 1:n if islin[i] == l && iseq[i] == e])
    end
    return order   # order[k] = natural index of reference row k
end

function check(name; tol = 1e-6)
    msgs = String[]
    bad(s) = push!(msgs, s)
    ref = CUTEstModel(name)
    try
        local m, perm
        try
            P = sif_decode(joinpath(ENV["MASTSIF"], name * ".SIF"))
            perm = cutest_perm(P)
            m = ExaModelsSIF.build_exa(P)
        catch err
            bad("DECODE ERROR: " * sprint(showerror, err))
            return msgs
        end
        ref.meta.nvar == m.meta.nvar ||
            bad("nvar: ref $(ref.meta.nvar) vs $(m.meta.nvar)")
        ref.meta.ncon == m.meta.ncon ||
            bad("ncon: ref $(ref.meta.ncon) vs $(m.meta.ncon)")
        isempty(msgs) || return msgs
        for (fld, a, b) in (("x0", ref.meta.x0, m.meta.x0),
                            ("lvar", ref.meta.lvar, m.meta.lvar),
                            ("uvar", ref.meta.uvar, m.meta.uvar),
                            ("lcon", ref.meta.lcon, m.meta.lcon[perm]),
                            ("ucon", ref.meta.ucon, m.meta.ucon[perm]))
            e = relerr(replace(a, Inf => 1e30, -Inf => -1e30),
                       replace(b, Inf => 1e30, -Inf => -1e30))
            e <= tol || bad("$fld mismatch (relerr $e)")
        end
        Random.seed!(7)
        pts = [("x0", copy(ref.meta.x0))]
        push!(pts, ("pert", clamp.(ref.meta.x0 .+ 0.01 .* randn(ref.meta.nvar),
                                   ref.meta.lvar, ref.meta.uvar)))
        for (tag, xp) in pts
            fa = obj(ref, xp); fb = obj(m, xp)
            e = abs(fa - fb) / (1.0 + max(abs(fa), abs(fb)))
            e <= tol || bad("obj@$tag: ref $fa vs $fb")
            e = relerr(grad(ref, xp), grad(m, xp))
            e <= tol || bad("grad@$tag (relerr $e)")
            if ref.meta.ncon > 0
                e = relerr(cons(ref, xp), cons(m, xp)[perm])
                e <= tol || bad("cons@$tag (relerr $e)")
                Ja = spjac(ref, xp); Jb = spjac(m, xp)[perm, :]
                e = norm(Ja - Jb, Inf) / (1.0 + max(norm(Ja, Inf), norm(Jb, Inf)))
                e <= tol || bad("jac@$tag (relerr $e)")
            end
            y = ref.meta.ncon > 0 ? randn(ref.meta.ncon) : Float64[]
            ym = zeros(length(y)); ym[perm] = y
            Ha = sphess(ref, xp, y); Hb = sphess(m, xp, ym)
            e = norm(Ha - Hb, Inf) / (1.0 + max(norm(Ha, Inf), norm(Hb, Inf)))
            e <= tol || bad("hess@$tag (relerr $e)")
        end
    catch err
        bad("ERROR: " * sprint(showerror, err))
    finally
        finalize(ref)
    end
    return msgs
end

function main(names = PILOT)
    npass = 0
    for name in names
        msgs = check(name)
        if isempty(msgs)
            npass += 1
            @printf("%-10s PASS\n", name)
        else
            @printf("%-10s FAIL  %s\n", name, msgs[1])
            for s in msgs[2:end]
                println(" "^16, s)
            end
        end
    end
    @printf("\n%d / %d pilot problems PASS\n", npass, length(names))
    return npass
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
