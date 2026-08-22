# Tests of utility functions

# Mainly just for "using"
include("common.jl")

const string_target = "Hello, JetReconstruction!"

struct IntegerOnly{T <: Integer} end

@testset "File streaming" begin
    for fname in ["file.txt", "file.txt.gz", "file.txt.zst"]
        fname = joinpath(@__DIR__, "data", fname)
        string = readline(JetReconstruction.open_with_stream(fname))
        @test string == string_target
    end
end

@testset "Return type concretization" begin
    @test JetReconstruction.concretize_return_type(PseudoJet, Float64) === PseudoJet
    @test JetReconstruction.concretize_return_type(Vector, Float32) === Vector{Float32}
    @test JetReconstruction.concretize_return_type(Vector{Float64}, Float32) ===
          Vector{Float64}

    @test_throws ArgumentError JetReconstruction.concretize_return_type(Array, Float64)
    @test_throws ArgumentError JetReconstruction.concretize_return_type(IntegerOnly,
                                                                        Float64)
    @test_throws ArgumentError JetReconstruction.concretize_return_type(Union{Int, Float64},
                                                                        Float64)
end
