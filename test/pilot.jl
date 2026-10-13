# test/pilot.jl — M1 pilot validation: compare every decoded model against the
# Fortran reference (CUTEst.jl) on metadata, values, and derivatives.
using ExaModelsSIF, ExaModels, CUTEst, NLPModels
using Random, Printf, LinearAlgebra, SparseArrays

ENV["MASTSIF"] = get(ENV, "MASTSIF", expanduser("~/git/pkg/cutest-sif"))

const PILOT = ["ROSENBR", "HS1", "HS2", "HS3", "HS4", "HS5", "HS6", "HS7",
    "HS10", "HS11", "HS14", "HS21", "HS26", "HS28", "HS35", "HS45",
    "BOOTH", "ZANGWIL3", "TAME", "ALSOTAME", "DEGDIAG", "GENHS28",
    # second wave: wide groups, variable/constraint reordering, var-SCALE,
    # conditional temporaries, Z-cards, first-set rule, param-operand codes
    "BROWNAL", "MEYER3", "DECONVC", "AVGASB", "CATENA", "GPP", "ARGLINA",
    "LUKVLE14"]

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

"Variable names as SIFDecode reports them (flattened array names)."
function refvarnames(ref)
    n = ref.meta.nvar
    status = zeros(Cint, 1)
    vname = Matrix{Cchar}(undef, 10, n)
    CUTEst.varnames(Float64, ref.libsif, status, Cint[n], vname)
    return [String(strip(join(Char(UInt8(c)) for c in vname[:, i]))) for i in 1:n]
end

"Constraint names as SIFDecode reports them (row order of the reference)."
function refconnames(ref)
    n = ref.meta.nvar; m = ref.meta.ncon
    status = zeros(Cint, 1)
    pname = Vector{Cchar}(undef, 10)
    vname = Matrix{Cchar}(undef, 10, n)
    cname = Matrix{Cchar}(undef, 10, m)
    CUTEst.cnames(Float64, ref.libsif, status, Cint[n], Cint[m], pname, vname, cname)
    return [String(strip(join(Char(UInt8(c)) for c in cname[:, i]))) for i in 1:m]
end

function check(name; tol = 1e-6)
    msgs = String[]
    bad(s) = push!(msgs, s)
    ref = CUTEstModel(name)
    try
        local m, perm, vp
        try
            P = sif_decode(joinpath(ENV["MASTSIF"], name * ".SIF"))
            m = ExaModelsSIF.build_exa(P)
            # CUTEst reorders both constraints and variables; map both by the
            # exact names SIFDecode reports
            connames = [gn for gn in P.grouporder if P.groups[gn].kind != 'N']
            rowof = Dict(gn => i for (i, gn) in enumerate(connames))
            perm = ref.meta.ncon > 0 ? [rowof[s] for s in refconnames(ref)] : Int[]
            vp = [P.varidx[s] for s in refvarnames(ref)]
        catch err
            bad("DECODE ERROR: " * sprint(showerror, err))
            return msgs
        end
        ref.meta.nvar == m.meta.nvar ||
            bad("nvar: ref $(ref.meta.nvar) vs $(m.meta.nvar)")
        ref.meta.ncon == m.meta.ncon ||
            bad("ncon: ref $(ref.meta.ncon) vs $(m.meta.ncon)")
        isempty(msgs) || return msgs
        for (fld, a, b) in (("x0", ref.meta.x0, m.meta.x0[vp]),
                            ("lvar", ref.meta.lvar, m.meta.lvar[vp]),
                            ("uvar", ref.meta.uvar, m.meta.uvar[vp]),
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
            xm = zeros(length(xp)); xm[vp] = xp
            fa = obj(ref, xp); fb = obj(m, xm)
            e = abs(fa - fb) / (1.0 + max(abs(fa), abs(fb)))
            e <= tol || bad("obj@$tag: ref $fa vs $fb")
            e = relerr(grad(ref, xp), grad(m, xm)[vp])
            e <= tol || bad("grad@$tag (relerr $e)")
            if ref.meta.ncon > 0
                e = relerr(cons(ref, xp), cons(m, xm)[perm])
                e <= tol || bad("cons@$tag (relerr $e)")
                Ja = spjac(ref, xp); Jb = spjac(m, xm)[perm, vp]
                e = norm(Ja - Jb, Inf) / (1.0 + max(norm(Ja, Inf), norm(Jb, Inf)))
                e <= tol || bad("jac@$tag (relerr $e)")
            end
            y = ref.meta.ncon > 0 ? randn(ref.meta.ncon) : Float64[]
            ym = zeros(length(y)); ym[perm] = y
            Ha = sphess(ref, xp, y)
            Hbf = sphess(m, xm, ym)
            Hbp = (Hbf + transpose(Hbf))[vp, vp]   # permute full, re-canonicalize
            Hap = Ha + transpose(Ha)
            e = norm(Hap - Hbp, Inf) / (1.0 + max(norm(Hap, Inf), norm(Hbp, Inf)))
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
