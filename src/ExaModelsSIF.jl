module ExaModelsSIF

using ExaModels

include("fortran.jl")
include("oracle_fallback.jl")
include("parser.jl")
include("builder.jl")

export sif_decode, sif_model

end # module
