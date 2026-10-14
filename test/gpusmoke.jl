# test/gpusmoke.jl — GPU backend smoke test: build each problem with
# CUDABackend and compare obj/grad/cons/jac/hess against the (validated)
# CPU model. Needs a CUDA-capable node and the env in the Pkg.activate line
# (develop ExaModels + ExaModelsSIF, add CUDA/KernelAbstractions/NLPModels).
# 13/13 PASS on 2x GV100 (2026-10-14); nonfinite patterns must match, finite
# entries to 1e-6 (grad/cons/jac) and 1e-5 x max|H| (hess, GPU sum order).
using Pkg; Pkg.activate(expanduser("~/git/scratch/gpuenv"))
using ExaModelsSIF, ExaModels, CUDA, NLPModels, SparseArrays, LinearAlgebra
const B = get(ENV, "MASTSIF", "/home/sushin/.julia/artifacts/1fa01d4c3b344b8cd32b8ef78c0abeac9bb9cce2")
# per-entry relative diff; CPU-nonfinite entries are compared by PATTERN
# (GPU must be nonfinite there too), matching the validation convention
function rel(a, b)
    any(isfinite.(a) .!= isfinite.(b)) && return Inf
    m = isfinite.(a)
    maximum(abs.(a[m] .- b[m]) ./ (1.0 .+ max.(abs.(a[m]), abs.(b[m]))); init = 0.0)
end
function smoke(name; params = Dict{String,Any}())
    P = ExaModelsSIF.sif_decode(joinpath(B, name * ".SIF"); params)
    mc = ExaModelsSIF.build_exa(P)
    mg = ExaModelsSIF.build_exa(P; backend = CUDABackend())
    x = copy(mc.meta.x0)
    xg = CuArray(x)
    msgs = String[]
    e = abs(obj(mc, x) - obj(mg, xg)) / (1 + abs(obj(mc, x)))
    e <= 1e-12 || push!(msgs, "obj $e")
    e = rel(grad(mc, x), Array(grad(mg, xg)))
    e <= 1e-6 || push!(msgs, "grad $e")
    if mc.meta.ncon > 0
        e = rel(cons(mc, x), Array(cons(mg, xg)))
        e <= 1e-6 || push!(msgs, "cons $e")
        rc, cc = jac_structure(mc); v = jac_coord(mc, x)
        Jc = sparse(rc, cc, v, mc.meta.ncon, mc.meta.nvar)
        rg = CuArray{Int}(undef, mg.meta.nnzj); cg = CuArray{Int}(undef, mg.meta.nnzj)
        jac_structure!(mg, rg, cg); vg = jac_coord(mg, xg)
        Jg = sparse(Array(rg), Array(cg), Array(vg), mc.meta.ncon, mc.meta.nvar)
        (nnz(Jc - Jg) == 0 || maximum(abs, Jc - Jg) <= 1e-6) || push!(msgs, "jac")
    end
    y = mc.meta.ncon > 0 ? ones(mc.meta.ncon) : Float64[]
    rh, ch = hess_structure(mc); vh = hess_coord(mc, x, y)
    Hc = sparse(rh, ch, vh, mc.meta.nvar, mc.meta.nvar)
    rhg = CuArray{Int}(undef, mg.meta.nnzh); chg = CuArray{Int}(undef, mg.meta.nnzh)
    hess_structure!(mg, rhg, chg)
    vhg = hess_coord(mg, xg, CuArray(y))
    Hg = sparse(Array(rhg), Array(chg), Array(vhg), mc.meta.nvar, mc.meta.nvar)
    d = Hc - Hg
    sc = 1.0 + maximum(abs, Hc; init = 0.0)
    (nnz(d) == 0 || maximum(abs, d) / sc <= 1e-5) || push!(msgs, "hess " * string(maximum(abs, d) / sc))
    println(rpad(name, 11), isempty(msgs) ? "GPU PASS" : "GPU FAIL: " * join(msgs, " | "))
end
for n in ("ROSENBR", "HS14", "DECONVC", "BROWNAL", "GULF", "MANCINO", "LHAIFAM",
          "DALLASM", "CHEBYQAD", "AVGASB", "LUKVLE14", "DEGENQP")
    try
        smoke(n)
    catch e
        println(rpad(n, 11), "GPU ERR: ", first(sprint(showerror, e), 90))
    end
end
smoke("FMINSURF"; params = Dict("P" => 16))
