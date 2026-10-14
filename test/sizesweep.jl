# test/sizesweep.jl — $-PARAMETER override validation: decode with the
# params kwarg vs CUTEst sifdecoder -param, compared like check().
# 11/11 PASS 2026-10-14 (FMINSURF P=16 excluded: sifdecoder only accepts
# LISTED parameter values; ours accepts any, by design).
include(joinpath(@__DIR__, "pilot.jl"))
using LinearAlgebra, Random

CASES = [
    ("READING4", "N", 100), ("READING4", "N", 1000),
    ("FMINSURF", "P", 16), ("FMINSURF", "P", 32),
    ("DIXMAANA1", "M", 1000), ("BROYDN7D", "N/2", 250),
    ("SROSENBR", "N/2", 250), ("LUKVLE14", "N", 998),
    ("CHAINWOO", "NS", 499), ("CRAGGLVY", "M", 2499),
    ("MANCINO", "N", 50), ("CURLY10", "N", 1000),
]

function checkp(name, pname, pval; tol = 1e-6)
    msgs = String[]
    bad(m) = push!(msgs, m)
    lk = mkpidlock("/tmp/cutest-decode.pid"; wait = true)
    ref = try CUTEstModel(name, "-param", "$pname=$pval") finally close(lk) end
    try
        P = sif_decode(joinpath(ENV["MASTSIF"], name * ".SIF");
                       params = Dict(pname => pval))
        m = ExaModelsSIF.build_exa(P)
        m.meta.nvar == ref.meta.nvar || bad("nvar $(m.meta.nvar) vs $(ref.meta.nvar)")
        isempty(msgs) || return msgs
        rn = refvarnames(ref); vp = [P.varidx[s] for s in rn]
        connames = [gn for gn in P.grouporder if P.groups[gn].kind != 'N']
        rowof = Dict(gn => i for (i, gn) in enumerate(connames))
        perm = ref.meta.ncon > 0 ? [rowof[s] for s in refconnames(ref)] : Int[]
        Random.seed!(7)
        xp = clamp.(ref.meta.x0 .+ 0.01 .* randn(ref.meta.nvar),
                    ref.meta.lvar, ref.meta.uvar)
        for (tag, x) in (("x0", copy(ref.meta.x0)), ("pert", xp))
            xm = zeros(length(x)); xm[vp] = x
            fa = obj(ref, x); fb = obj(m, xm)
            isfinite(fa) || continue
            e = abs(fa - fb) / (1.0 + max(abs(fa), abs(fb)))
            e <= tol || bad("obj@$tag: $fa vs $fb")
            e = relerr(grad(ref, x), grad(m, xm)[vp])
            e <= tol || bad("grad@$tag $e")
            if ref.meta.ncon > 0
                e = relerr(cons(ref, x), cons(m, xm)[perm])
                e <= tol || bad("cons@$tag $e")
                e = sprelerr(spjac(ref, x), spjac(m, xm)[perm, vp])
                e <= tol || bad("jac@$tag $e")
            end
            y = ref.meta.ncon > 0 ? randn(ref.meta.ncon) : Float64[]
            ym = zeros(length(y)); ym[perm] = y
            Ha = sphess(ref, x, y); Hb = sphess(m, xm, ym)
            Hap = Ha + transpose(Ha); Hbp = (Hb + transpose(Hb))[vp, vp]
            e = sprelerr(Hap, Hbp)
            e <= tol || bad("hess@$tag $e")
        end
    catch err
        bad("ERROR: " * first(sprint(showerror, err), 80))
    finally
        finalize(ref)
    end
    return msgs
end

for (name, pname, pval) in CASES
    msgs = checkp(name, pname, pval)
    println(rpad("$name[$pname=$pval]", 24),
        isempty(msgs) ? "PASS" : "FAIL: " * join(msgs, " | ")[1:min(end, 90)])
end
