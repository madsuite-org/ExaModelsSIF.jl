# oracle_fallback.jl — external-Fortran tier: problems whose element/group
# functions are hand-coded Fortran appended to the SIF cannot be decoded into
# expressions; interface the compiled Fortran instead (Sungho, 2026-10-13).
# CUTEst.jl compiles exactly that Fortran; its callbacks are wrapped in
# ExaModels oracles so sif_model covers these problems through one entry point
# (CPU callbacks; no SIMD/GPU for this tier).

import CUTEst
import NLPModels

"Raised at decode when a type's function body calls external Fortran."
struct ExternalFortranError <: Exception
    tname::String
    fn::Symbol
end
Base.showerror(io::IO, e::ExternalFortranError) = print(io,
    "type $(e.tname) uses external Fortran function $(e.fn) (oracle tier)")

"""
    cutest_oracle_model(path; T = Float64) -> ExaModel

Build an ExaModel whose objective and constraints are ExaModels oracles backed
by a `CUTEstModel` (which compiles the SIF's embedded Fortran). Variable and
constraint order are CUTEst's own.
"""
function cutest_oracle_model(path::AbstractString; T = Float64)
    dir = dirname(abspath(path))
    name = replace(basename(path), r"\.SIF$" => "")
    ENV["MASTSIF"] = dir
    ref = CUTEst.CUTEstModel(name)
    meta = ref.meta
    n, m = meta.nvar, meta.ncon

    c = ExaModels.ExaCore(T)
    c, x = ExaModels.add_var(c, n; start = meta.x0, lvar = meta.lvar, uvar = meta.uvar)

    hrows, hcols = NLPModels.hess_structure(ref)
    y0 = zeros(m)
    sor = ExaModels.ScalarNonlinearOracle(;
        nvar = n,
        f = x_ -> (ref_ = ref; NLPModels.obj(ref_, x_)),
        grad! = (g, x_) -> NLPModels.grad!(ref, x_, g),
        hvp! = (Hv, x_, v) -> NLPModels.hprod!(ref, x_, y0, v, Hv; obj_weight = 1.0),
        hess_rows = copy(hrows), hess_cols = copy(hcols),
        adapt = Val(true),
    )
    c = ExaModels.objective(c, sor)

    if m > 0
        jrows, jcols = NLPModels.jac_structure(ref)
        vor = ExaModels.VectorNonlinearOracle(;
            nvar = n, ncon = m,
            jac_rows = copy(jrows), jac_cols = copy(jcols),
            hess_rows = copy(hrows), hess_cols = copy(hcols),
            lcon = meta.lcon, ucon = meta.ucon,
            f! = (cv, x_) -> NLPModels.cons!(ref, x_, cv),
            jac! = (vals, x_) -> NLPModels.jac_coord!(ref, x_, vals),
            hess! = (vals, x_, y) -> NLPModels.hess_coord!(ref, x_, y, vals; obj_weight = 0.0),
            adapt = Val(true),
        )
        c = ExaModels.constraint(c, vor)
    end
    return ExaModels.ExaModel(c)
end
