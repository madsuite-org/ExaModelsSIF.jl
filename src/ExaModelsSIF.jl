module ExaModelsSIF

using ExaModels

include("fortran.jl")
include("parser.jl")
include("builder.jl")

export sif_decode, sif_model, ExternalFortranError

end # module
