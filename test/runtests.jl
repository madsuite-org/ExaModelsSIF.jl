using Test
include("eltype.jl")
include("pilot.jl")
@test main() == length(PILOT)
