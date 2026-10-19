# ExaModelsSIF.jl

[![CI](https://github.com/madsuite-org/ExaModelsSIF.jl/actions/workflows/test.yml/badge.svg?branch=master)](https://github.com/madsuite-org/ExaModelsSIF.jl/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/madsuite-org/ExaModelsSIF.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/madsuite-org/ExaModelsSIF.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

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

The decoder is implemented from the SIF reference report [3], which
SIFDecode distributes.

[S2MPJ](https://github.com/GrattonToint/S2MPJ) [5] also translates the SIF
collection into native code, producing problem files for Matlab, Python and
Julia that evaluate values and derivatives without Fortran. ExaModelsSIF
targets a different form: it emits ExaModels' SIMD representation [6, 7], in
which instances of one element or group type share one compiled kernel, so a
problem evaluates on a GPU and derivatives come from ExaModels' automatic
differentiation.

## Coverage and validation

The CUTEst collection distributed with CUTEst.jl [8] v1.4.0 has 1548
problems. 1514 of them decode and build, and each is compared against the
compiled Fortran [1] through CUTEst.jl: objective, gradient, constraints, Jacobian and
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

## References

1. N. I. M. Gould, D. Orban, Ph. L. Toint, *CUTEst: a Constrained and Unconstrained Testing Environment with safe threads for mathematical optimization*, Computational Optimization and Applications 60(3), 545-557, 2015. [doi:10.1007/s10589-014-9687-3](https://doi.org/10.1007/s10589-014-9687-3)
2. A. R. Conn, N. I. M. Gould, Ph. L. Toint, *LANCELOT: A Fortran Package for Large-Scale Nonlinear Optimization (Release A)*, Springer Series in Computational Mathematics, Springer, 1992. [doi:10.1007/978-3-662-12211-2](https://doi.org/10.1007/978-3-662-12211-2)
3. A. R. Conn, N. I. M. Gould, D. Orban, Ph. L. Toint, *The SIF Reference Report (revised version)*, distributed with SIFDecode. [sif.pdf](https://github.com/ralna/SIFDecode/blob/master/doc/pdf/sif.pdf)
4. SIFDecode, the SIF decoder of CUTEst. [github.com/ralna/SIFDecode](https://github.com/ralna/SIFDecode)
5. S. Gratton, Ph. L. Toint, *S2MPJ and CUTEst optimization problems for Matlab, Python and Julia*, Optimization Methods and Software 40(4), 871-903, 2025. [doi:10.1080/10556788.2025.2490640](https://doi.org/10.1080/10556788.2025.2490640), [arXiv:2407.07812](https://arxiv.org/abs/2407.07812)
6. S. Shin, M. Schanen, F. Pacaud, A. Montoison, M. Anitescu, *ExaModels.jl: An Algebraic Modeling System for Nonlinear Programming on GPUs*, 2026. [arXiv:2608.16265](https://arxiv.org/abs/2608.16265)
7. S. Shin, M. Anitescu, F. Pacaud, *Accelerating optimal power flow with GPUs: SIMD abstraction of nonlinear programs and condensed-space interior-point methods*, Electric Power Systems Research 236, 110651, 2024. [doi:10.1016/j.epsr.2024.110651](https://doi.org/10.1016/j.epsr.2024.110651)
8. A. Montoison, D. Orban, A. Soares Siqueira, and contributors, *CUTEst.jl: Julia's CUTEst interface*. [doi:10.5281/zenodo.1188851](https://doi.org/10.5281/zenodo.1188851)
9. D. Orban, A. Soares Siqueira, and contributors, *NLPModels.jl: Data Structures for Optimization Models*. [doi:10.5281/zenodo.2558627](https://doi.org/10.5281/zenodo.2558627)

## Acknowledgements

The SIF format comes from LANCELOT [2], and the problems are those of the
CUTEst collection [1], maintained by the CUTEst developers. SIFDecode [4],
compiled and called through CUTEst.jl [8], is the reference every model here
is validated against. Models implement the NLPModels.jl [9] interface through
ExaModels.

## License

MIT; see [LICENSE](LICENSE). No CUTEst Fortran or C source is included or
distributed. The SIF problem files are distributed by the CUTEst project under
their own terms.
