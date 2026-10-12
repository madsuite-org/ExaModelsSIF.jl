using Test
include("pilot.jl")
@test main() == length(PILOT)
