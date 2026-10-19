# ExaModelsSIF.jl

ExaModelsSIF decodes optimization problems written in SIF, the Standard Input
Format of the [CUTEst](https://github.com/ralna/CUTEst) test set, directly into
[ExaModels](https://github.com/exanauts/ExaModels.jl) models. The resulting
models are evaluated and differentiated in Julia, on the CPU or on any GPU
backend ExaModels supports, with no Fortran compiler and no SIF preprocessor.

```julia
using ExaModelsSIF

m = sif_model("ROSENBR.SIF")                      # an ExaModels.ExaModel
m = sif_model("ROSENBR.SIF"; T = Float32)         # any floating-point type
m = sif_model("DIXMAANA1.SIF"; backend = CUDABackend())   # on a GPU (using CUDA)

P = sif_decode("DIXMAANA1.SIF"; params = Dict("M" => 1000))  # the decoded problem
```

For the CUTEst collection itself, with problem selection by name and
classification, use [ExaCUTEst.jl](https://github.com/madsuite-org/ExaCUTEst.jl),
which is built on this package.

## How it works

SIF describes a problem as group partially separable: each objective term and
each constraint is a group function applied to a scaled sum of linear terms and
nonlinear elements, and all elements or groups of one type share one algebraic
expression. The decoder

1. interprets the SIF data sections (fixed-field cards, `DO` loops, integer and
   real parameter arithmetic, array names) into a fully instantiated
   representation of the problem;
2. translates the Fortran expressions of the element and group sections into
   Julia expressions; derivatives come from ExaModels' automatic
   differentiation, so the hand-coded `G` and `H` lines are not used;
3. groups elements and groups by type signature and emits one ExaModels
   generator per group of instances, which is the SIMD form ExaModels compiles
   into a single kernel.

Constraint rows keep SIF group order. The group constant is folded into the
constraint value, following the CUTEst convention. Internal variables are
substituted as linear combinations, and the `QUADRATIC` section is emitted as
Hessian terms. Groups too wide to trace as one kernel are lowered through
buffered subexpressions and a reduction.

The decoder is implemented from the published
[SIF reference](https://www.numerical.rl.ac.uk/lancelot/sif/).

## Coverage and validation

The CUTEst collection distributed with CUTEst.jl v1.4.0 has 1548 problems.
1514 of them decode and build, and each is compared against the compiled
Fortran through CUTEst.jl: objective, gradient, constraints, Jacobian and
Lagrangian Hessian at the starting point and at a perturbed point, together
with sizes, starting point and bounds. All 1514 agree.

The other 34 append hand-written Fortran routines to their SIF files, which a
decoder cannot translate. For these, `sif_decode` raises `ExternalFortranError`
by default, and `external = :defer` returns the decoded structure with those
element types left for the caller to supply. ExaCUTEst.jl provides native Julia
implementations of all 34.

A few disagreements with the reference turned out to be defects in the
reference problems, reported upstream as
[SIFDecode#55](https://github.com/ralna/SIFDecode/issues/55) and
[#56](https://github.com/ralna/SIFDecode/issues/56). In several others the
SIF file's hand-coded Hessian contradicts a finite difference of its own
gradient; the validation adjudicates those against that finite difference.

`test/pilot.jl` runs the comparison against CUTEst.jl, and `test/eltype.jl`
checks that a model built at element type `T` computes in `T` throughout.

## Installation

ExaModelsSIF is not yet registered, and it needs the `sif-integration`
branch of [madsuite-org/ExaModels.jl](https://github.com/madsuite-org/ExaModels.jl),
which carries changes not yet in a released ExaModels. Install that branch
first, since Pkg does not follow a dependency's own `[sources]`:

```julia
using Pkg
Pkg.add(url = "https://github.com/madsuite-org/ExaModels.jl", rev = "sif-integration")
Pkg.add(url = "https://github.com/madsuite-org/ExaModelsSIF.jl")
```

No CUTEst Fortran or C source is included or distributed.
