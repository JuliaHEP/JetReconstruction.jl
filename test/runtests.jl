using JetReconstruction
using Test

"""
A function for test comparison
"""
function arrcompare(y, yt; eps=0.01)
    for i in 1:length(y)
        if !(sum(abs.(y[i] .- yt[i]) .< eps) == length(y[i]))
            return false
        end
    end
    true
end

const NUMBER_OF_TESTS = 11 # number of test files in the data folder

@testset "JetReconstruction.jl" begin
    # @test anti_kt(X) == Y

    # additional simple test
    @test arrcompare(
        anti_kt([
            [100, 99.0, 0.1, 0],
            [5.0, 4.0, -0.1, 0],
            [99., -99, 0., 0.0]
        ], R=0.7)[1],

        [[105.0, 103.0, 0.0, 0.0], [99.0, -99.0, 0.0, 0.0]]
    )

    # actual testing
    for i in 1:NUMBER_OF_TESTS
        istr = string(i)

        @test arrcompare(
                sort(anti_kt(loadjets(
                    "data/"*istr*".dat", constructor=(x,y,z,E)->[E,x,y,z]
                ), R=1)[1], lt=(a,b) -> (pt(a) > pt(b))),

                loadjets(
                    "data/"*istr*"-fj-result.dat", constructor=(x,y,z,E)->[E,x,y,z]
                )
            )
    end
end
