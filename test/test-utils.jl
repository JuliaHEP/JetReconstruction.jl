# Tests of utility functions

# Mainly just for "using"
include("common.jl")

const string_target = "Hello, JetReconstruction!"

@testset "File streaming" begin
    for fname in ["file.txt", "file.txt.gz", "file.txt.zst"]
        fname = joinpath(@__DIR__, "data", fname)
        string = readline(JetReconstruction.open_with_stream(fname))
        @test string == string_target
    end
end
