using JetReconstruction
using Test

@testset "JetReconstruction.jl" begin
    # @test anti_kt(X) == Y
    @test anti_kt([
        [100, 99.0, 0.1, 0],
        [5.0, 4.0, -0.1, 0],
        [99., -99, 0., 0.0]
    ], R=0.7)[1] == [[105.0, 103.0, 0.0, 0.0], [99.0, -99.0, 0.0, 0.0]]

    #=for i in 1:10
        istr = string(i)

        @test sort(anti_kt(
                loadjets(
                    "data/"*istr*".dat", constructor=(x,y,z,E)->[E,x,y,z]
                ), R=1)[1], lt=(a,b) -> (pt(a) > pt(b))) == loadjets(
                    "data/"*istr*"-fj-result.dat", constructor=(x,y,z,E)->[E,x,y,z]
                )
    end=#
end
