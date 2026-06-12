# Tests of utility functions

include("common.jl")

const string_target = "Hello, JetReconstruction!"

@testset "File streaming" begin
    for fname in ["file.txt", "file.txt.gz", "file.txt.zst"]
        fname = joinpath(@__DIR__, "data", fname)
        string = readline(JetReconstruction.open_with_stream(fname))
        @test string == string_target
    end
end

@testset "Return type concretization" begin
    @test JetReconstruction.concretize_return_type(EEJet, Float64) == EEJet{Float64}
    @test_throws ArgumentError _=JetReconstruction.concretize_return_type(Array, Float64)
end
