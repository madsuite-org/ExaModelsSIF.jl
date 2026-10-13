# ExaModelsSIF.jl

A pure-Julia SIF (Standard Input Format) decoder that converts CUTEst test
problems directly into [ExaModels](https://github.com/exanauts/ExaModels.jl)
models, so the CUTEst collection can be evaluated and differentiated natively
on CPU and GPU backends without Fortran.

```julia
using ExaModelsSIF
m = sif_model("path/to/ROSENBR.SIF")           # -> ExaModels.ExaModel
P = sif_decode("path/to/ROSENBR.SIF")          # -> instantiated problem IR
```

## How it works

SIF is group-partially-separable: every function is
`GF(scale^-1 * (sum_j a_j x_j + sum_e w_e f_e(x_e) - b))`, and all
elements/groups of one TYPE share one algebraic expression. The decoder:

1. interprets the SIF data sections (fixed-field cards, DO loops, integer/real
   parameter arithmetic, array names) into a fully instantiated IR;
2. translates the restricted-Fortran `F =` expressions of the
   ELEMENTS/GROUPS sections into Julia ASTs (derivatives come from ExaModels'
   own AD, so `G`/`H` lines are ignored);
3. buckets elements and groups by type signature and emits one ExaModels
   generator per bucket — the SIMD abstraction — via `add_var`/`add_obj`/
   `add_con`/`add_con!`.

Constraint rows keep natural SIF group order. The group constant is folded
into the constraint value (CUTEst convention), internal variables (`R` range
transformations) are substituted as linear combinations at trace time, and
the `QUADRATIC` section is emitted as Hessian-entry terms.

Clean-room implementation from the published SIF reference
(https://www.numerical.rl.ac.uk/lancelot/sif/); no SIFDecode/CUTEst source
was consulted.

## Validation

`test/pilot.jl` compares every decoded model against CUTEst.jl (Fortran
reference): nvar/ncon, x0, variable and constraint bounds, objective,
gradient, constraints, Jacobian, and Lagrangian Hessian at x0 and a perturbed
point. 22/22 pilot problems pass. `test/sweep.jl` decodes + builds the whole
local collection and writes `sweep_results.tsv`.

## Reproducibility

All scripts run under the project environment (`julia --project=.`); the
committed `Manifest.toml` pins every dependency. `ExaModels` is consumed as
a path source (`~/git/pkg/ExaModels.jl`, shared across the cluster via NFS);
the validated revision is madsuite-org/ExaModels.jl branch compile-time-prototypes, commit `fe407fc0968ac52762965f57a1f10767ad9da5d9` (worktree ~/git/worktrees/ExaModels-ctp). See STATUS.md.
