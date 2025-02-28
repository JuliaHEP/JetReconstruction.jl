# Tests of algorithm/power consistency checks

include("common.jl")

@testset "Algorithm/power consistency" begin
    @test JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.AntiKt,
                                                              p = -1)
    @test JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.CA,
                                                              p = 0)
    @test JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.Kt,
                                                              p = 1)

    @test JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.AntiKt,
                                                              p = nothing)
    @test JetReconstruction.check_algorithm_power_consistency(algorithm = nothing,
                                                              p = -1)

    @test_throws ArgumentError JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.AntiKt,
                                                                                   p = 0)
    @test_throws ArgumentError JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.Kt,
                                                                                   p = 1.5)

    @test JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.GenKt,
                                                              p = 1.5)
    @test JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.GenKt,
                                                              p = -0.5)

    @test_throws ArgumentError JetReconstruction.check_algorithm_power_consistency(algorithm = JetAlgorithm.GenKt,
                                                                                   p = nothing)
end
