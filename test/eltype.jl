# A decoded model must be pure in its element type: every number it returns,
# and every number it computes with, carries `T`.
#
# The visible symptom of a leak is `obj` returning a `Float64` from a
# `T = Float32` model, and it is a weak symptom in a telling way: `grad` and
# `cons` stay Float32 even when the leak is present, because they are written
# into the model's own typed arrays. So the array element types agree while
# the arithmetic that produced them ran in Float64 and was rounded once on the
# way out -- which is wrong in the direction that looks right, and on a GPU
# backend is Float64 arithmetic inside a Float32 kernel. `obj` is the only
# return value that is a bare scalar, so it is the one reading that exposes
# the data tables' type.
using ExaModelsSIF, ExaModels, NLPModels, Test
import CUTEst   # sets ENV["MASTSIF"] to its SIF collection unless already set

const ART = ENV["MASTSIF"]

# Chosen to cover the distinct emit paths in `build_exa`, since each builds its
# own data rows: a plain objective group, linear and nonlinear constraint rows,
# an objective whose groups carry a nontrivial group function, a wide group
# lowered through `add_reduce`, a quadratic section, and a problem with
# element parameters and injected temporaries.
const PROBLEMS = ("ROSENBR", "HS14", "HS35", "WATSON", "BROWNAL", "GULF",
                  "LUKVLE2", "CLEUVEN7")

@testset "element type purity" begin
    for nm in PROBLEMS
        path = joinpath(ART, nm * ".SIF")
        isfile(path) || continue
        @testset "$nm" begin
            for T in (Float32, Float64)
                m = ExaModelsSIF.sif_model(path; T)
                x = m.meta.x0
                @test eltype(x) === T
                @test typeof(obj(m, x)) === T
                @test eltype(grad(m, x)) === T
                @test eltype(m.meta.lvar) === T
                if m.meta.ncon > 0
                    @test eltype(cons(m, x)) === T
                    @test eltype(m.meta.lcon) === T
                end
            end
        end
    end
    # The Float64 case above is the control: it reads as pure whether or not
    # the conversion is in place, so on its own it cannot witness the property.
    @test length(PROBLEMS) == 8
end
