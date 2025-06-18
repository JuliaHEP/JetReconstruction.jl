# Tests of algorithm/power consistency checks

include("common.jl")

@testset "Algorithm/power consistency" begin
    @test JetReconstruction.get_algorithm_power(algorithm = JetAlgorithm.AntiKt,
                                                p = nothing) == -1
    @test JetReconstruction.get_algorithm_power(algorithm = JetAlgorithm.CA,
                                                p = nothing) == 0
    @test JetReconstruction.get_algorithm_power(algorithm = JetAlgorithm.Kt,
                                                p = nothing) == 1
    @test JetReconstruction.get_algorithm_power(algorithm = JetAlgorithm.Durham,
                                                p = nothing) == 1

    @test JetReconstruction.get_algorithm_power(algorithm = JetAlgorithm.GenKt,
                                                p = 1.5) == 1.5
    @test JetReconstruction.get_algorithm_power(algorithm = JetAlgorithm.GenKt,
                                                p = -0.5) == -0.5
    @test JetReconstruction.get_algorithm_power(algorithm = JetAlgorithm.EEKt,
                                                p = 1.5) == 1.5
    @test JetReconstruction.get_algorithm_power(algorithm = JetAlgorithm.EEKt,
                                                p = -0.5) == -0.5

    @test_throws ArgumentError JetReconstruction.get_algorithm_power(algorithm = JetAlgorithm.GenKt,
                                                                     p = nothing)
    @test_throws ArgumentError JetReconstruction.get_algorithm_power(algorithm = JetAlgorithm.EEKt,
                                                                     p = nothing)
end
