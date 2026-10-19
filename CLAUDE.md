# ExaModelsSIF.jl developer guide

ExaModelsSIF decodes SIF files into ExaModels models: SIF text in, `ExaModel`
out. It has no dependency on CUTEst.jl, no Fortran toolchain requirement, and
no package extension for either. A problem whose SIF file appends external
Fortran raises `ExternalFortranError`; ExaCUTEst.jl provides native
implementations of those problems. CUTEst.jl is a test dependency only, as the
reference the decoder is validated against.

## Layout

```
src/
  ExaModelsSIF.jl   module; exports sif_decode, sif_model, ExternalFortranError
  parser.jl         fixed-field card tokenizer and data-section interpreter
                    (DO loops, parameter arithmetic, array names) producing a
                    fully instantiated problem (SIFProblem); the ELEMENTS and
                    GROUPS function-part parser; external-Fortran detection
  fortran.jl        Fortran expression language to Julia ASTs, evaluated at
                    ExaModels trace time; registered quiet intrinsics
  builder.jl        SIFProblem to ExaModel (build_exa, sif_model)
test/
  runtests.jl       eltype.jl, then pilot.jl over its problem list
  eltype.jl         a model built at element type T computes in T throughout
  pilot.jl          comparison against CUTEst.jl; check() and the
                    adjudication sets for defects in the reference
  sizesweep.jl      $-PARAMETER overrides against sifdecoder -param
  valsweep.jl       pilot.jl's comparison over the whole collection
  sweep.jl          decode and build over the whole collection, no reference
  gpusmoke.jl       GPU models against CPU models (needs CUDA; not in Pkg.test)
```

The decoder is implemented from the published SIF reference
(https://www.numerical.rl.ac.uk/lancelot/sif/). Derivatives come from
ExaModels' automatic differentiation; the hand-coded `G` and `H` lines of a SIF
file are not used.

## How a problem becomes an ExaModel

Elements and groups are bucketed by type signature, and each bucket becomes one
ExaModels generator over a concrete vector of named tuples, so instances of one
type share one compiled kernel. Nontrivial groups wider than `WIDE_LIN` linear
terms, or with more than `ELEM_BUDGET` element fields, lower their wide part
through buffered subexpressions and a reduction, so no traced kernel grows with
the width of a group. A buffered stage spans all instances of one element type;
stages live in the model's type, so their count must never scale with the
problem size.

Constraint rows keep SIF group order and variables keep SIF declaration order.
CUTEst reorders both (nonlinear variables first; equality and linear rows
first), so a comparison against CUTEst.jl permutes by SIF name. `pilot.jl`
reads the reference's names through `CUTEst.varnames`/`CUTEst.cnames`; never
guess the order.

## SIF semantics

These are measured against the Fortran reference; each was confirmed by the
named problem agreeing with it.

- Default variable bounds are `[0, Inf)`. Bound magnitudes `>= 1e20` are
  infinite, for variables and constraints.
- An `MI`/`XM`/`ZM` bound sets the lower bound to `-Inf`, and an upper bound
  still at the default `+Inf` drops to 0 (MATRIX2).
- A `BOUNDS 'DEFAULT'` card also governs variables created after the BOUNDS
  section, for example by `ZV` in ELEMENT USES (NUFFIELD).
- A `RANGES 'DEFAULT'` card leaves equality rows untouched and ranges the
  `G`/`L` rows (DEGENQP, DEGENQPC).
- Only the first named set applies in CONSTANTS, RANGES, BOUNDS and START
  POINT.
- The group constant is folded into the constraint value with a minus sign,
  and bounds are built around 0.
- A group `SCALE` divides the group value. Variable `SCALE` factors are
  LANCELOT hints and are ignored, since CUTEst evaluates the unscaled problem
  (MEYER3).
- Repeated cards naming the same group add their linear terms to it.
- `INTEGER`/`CONTINUOUS` markers are ignored: the problem is its continuous
  relaxation.
- The QUADRATIC section uses the Hessian-entry reading: an entry `(i, j, q)`
  contributes `q*x[i]*x[j]` off the diagonal and `q*x[i]^2/2` on it.
- A DO-loop variable keeps its last value after the loop.
- `**` with an integer-declared (`I`) temporary exponent is Fortran integer
  power, defined for a negative base, lowered to `_npow`; integer temporary
  assignments truncate (MANCINO).
- `LOG`, `DLOG` and `DSQRT` are quiet: they return NaN outside their domain
  rather than throwing, through registered functions (`_flog`, `_fsqrt`).
- An undecorated real literal in a function section is default `REAL` in the
  Fortran SIFDecode generates, so it is rounded through `Float32`;
  `D`-exponent literals stay double and integers stay integers (READING4,
  READING5).

## Disagreements with the reference

SIFDecode writes every double into `OUTSDIF.d` with 11 significant digits, so
the reference evaluates on truncated data while a decoded model uses the SIF's
values. This is the usual cause of a small, reproducible, structurally clean
disagreement. Its size depends on what the SIF does with the value: about
1e-11 passed straight through, up to about 1e-8 after a cancelling
subtraction, or a flat relative offset on every value a group scales. The
`x0` a `CUTEstModel` reports comes from the same file, so comparisons evaluate
both models at the reference's `x0`.

Five problems in the collection supply first derivatives only (a `1` in the
fourth classification field): BLEACHNG, CAR2, ORBIT2, SARO and SAROMM. All five
append external Fortran, so no decoded model is affected. Their reference
Hessians are absent or partial; adjudicate them against a central difference
of the reference's Lagrangian gradient.

`test/pilot.jl` names every problem whose reference is known to be wrong, with
the evidence beside each set:

- `HESS_REF_BUGS`: hand-coded `H` lines contradict a finite difference of the
  reference's own gradient while ours agrees;
- `DERIV_REF_BUGS`: reference `F` and `G` are inconsistent by construction
  (EQC);
- `VALUE_REF_BUGS`: the hand-coded value contradicts the SIF text (HIMMELP*);
- `REF_DECODE_BUGS`: SIFDecode decodes the file wrongly (LRW8A);
- `REF_NOISE`, `REF_NOISE_TOL`: the reference's own floating-point noise,
  adjudicated in `BigFloat`.

Adding a problem to one of these sets needs the same kind of evidence, usually
a `BigFloat` evaluation or a finite difference of the reference. Never loosen
a tolerance to cover a disagreement.

## Running the tests

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

Needs `gfortran`, since the reference builds `CUTEstModel`s. `Pkg.test()` uses
`test/Project.toml`; the scripts run in that environment directly:

```
julia --project=test test/pilot.jl
julia --project=test test/sizesweep.jl
bash test/valsweep_driver.sh 12     # whole collection, sharded and restartable
bash test/sweep_driver.sh 8         # decode and build only
```

Loading CUTEst.jl sets `ENV["MASTSIF"]` to the collection it ships unless it is
already set, so the decoder and the reference read the same files. The sweeps
read `SIF_DIR` or `MASTSIF` and write their results to the repository root,
where `.gitignore` covers them.

CUTEst.jl decodes and compiles every problem into one shared directory, so
concurrent `CUTEstModel` builds on one machine corrupt each other. Every build
in the test scripts holds `mkpidlock("/tmp/cutest-decode.pid"; wait = true)`;
keep that for any new one. A `CUTEstModel` must be finalized before another
with the same name is built.

`test/gpusmoke.jl` needs a CUDA device and an environment with CUDA; it is not
part of `Pkg.test()`. GPU models take device structure buffers
(`jac_structure!`/`hess_structure!` with device integer arrays), which is
ExaModels' calling convention.

## Dependencies

ExaModelsSIF needs the `sif-integration` branch of madsuite-org/ExaModels.jl,
declared in `[sources]` of both `Project.toml` and `test/Project.toml`. The
tracked `Manifest.toml` pins the revision tested.
