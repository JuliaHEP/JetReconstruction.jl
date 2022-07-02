# A temporary file to run some quick tests

using StaticArrays

include("src/JetReconstruction.jl")
using .JetReconstruction
# using JetReconstruction

# unrealistic tests
anti_kt([
    SVector(π, 0, 0, 1)
])

anti_kt([
    SVector(π, 0, 0, 1),
    SVector(π, 0, 2, 0),
    SVector(1, 1, 0, 0),
    SVector(1, 0, 1, 0),
])
