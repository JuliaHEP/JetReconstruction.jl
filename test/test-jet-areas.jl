using Test
using JetReconstruction

@testset "active jet areas with explicit ghosts" begin
    particles = [
        PseudoJet(100.0, 0.0, 0.0, 100.0; cluster_hist_index = 1)
    ]

    area_spec = GhostedAreaSpec(ghost_maxrap = 1.0,
                                ghost_area = 0.02,
                                seed = 12345)

    csa = jet_reconstruct_area(particles;
                               algorithm = JetAlgorithm.AntiKt,
                               R = 0.4,
                               area_spec = area_spec,
                               strategy = RecoStrategy.N2Plain,)

    jets = inclusive_jets(csa, PseudoJet; ptmin = 1.0)
    default_jets = inclusive_jets(csa; ptmin = 1.0)

    @test length(jets) == 1
    @test eltype(default_jets) == PseudoJet
    @test area(default_jets[1], csa) ≈ area(jets[1], csa)
    @test !is_pure_ghost(jets[1], csa)
    @test area(jets[1], csa) > 0.0
    @test area(jets[1], csa)≈π * 0.4^2 rtol=0.25
    @test total_area(csa) > 0.0
    @test has_explicit_ghosts(csa)
end

@testset "area metadata helpers and return conversions" begin
    particles = [
        PseudoJet(100.0, 0.0, 0.0, 100.0; cluster_hist_index = 1)
    ]

    csa = jet_reconstruct_area(particles;
                               algorithm = JetAlgorithm.AntiKt,
                               R = 0.4,
                               area_spec = GhostedAreaSpec(ghost_maxrap = 1.0,
                                                           ghost_area = 0.02,
                                                           seed = 12345),)

    jets = inclusive_jets(csa; ptmin = 1.0)
    jet = only(jets)

    @test area_4vector(jet, csa) isa PseudoJet
    @test area_error(jet, csa) == 0.0
    @test max_ghost_pt2(csa) > 0.0
    @test !has_dangerous_particles(csa)
    @test n_hard_particles(csa) == length(particles)

    lv_jets = inclusive_jets(csa, JetReconstruction.LorentzVector; ptmin = 1.0)
    @test lv_jets isa Vector
    @test eltype(lv_jets) <: JetReconstruction.LorentzVector

    lvc_jets = inclusive_jets(csa, JetReconstruction.LorentzVectorCyl; ptmin = 1.0)
    @test lvc_jets isa Vector
    @test eltype(lvc_jets) <: JetReconstruction.LorentzVectorCyl
end

@testset "pure ghost jets are filtered by default" begin
    particles = [
        PseudoJet(100.0, 0.0, 0.0, 100.0; cluster_hist_index = 1)
    ]

    area_spec = GhostedAreaSpec(ghost_maxrap = 1.0,
                                ghost_area = 0.05,
                                seed = 12345)

    csa = jet_reconstruct_area(particles;
                               algorithm = JetAlgorithm.AntiKt,
                               R = 0.4,
                               area_spec = area_spec,
                               strategy = RecoStrategy.N2Plain,)

    physical_jets = inclusive_jets(csa, PseudoJet; ptmin = 0.0)
    all_jets = inclusive_jets(csa, PseudoJet; ptmin = 0.0, include_pure_ghosts = true)

    @test length(all_jets) >= length(physical_jets)
    @test all(!is_pure_ghost(jet, csa) for jet in physical_jets)
end

@testset "ghost area validation and dangerous-particle metadata" begin
    particles = [
        PseudoJet(100.0, 0.0, 0.0, 100.0; cluster_hist_index = 1)
    ]

    @test_throws ArgumentError jet_reconstruct_area(particles;
                                                    algorithm = JetAlgorithm.AntiKt,
                                                    R = 0.4,
                                                    area_spec = GhostedAreaSpec(repeat = 2),)

    csa = jet_reconstruct_area(particles;
                               algorithm = JetAlgorithm.AntiKt,
                               R = 0.4,
                               area_spec = GhostedAreaSpec(ghost_maxrap = 1.0,
                                                           ghost_area = 0.05,
                                                           mean_ghost_pt = 1.0e-5,
                                                           seed = 12345),
                               preprocess = nothing,)

    @test has_dangerous_particles(csa)
end

@testset "area propagation follows clustering history" begin
    particles = [
        PseudoJet(100.0, 0.0, 0.0, 100.0; cluster_hist_index = 1),
        PseudoJet(-80.0, 0.0, 0.0, 80.0; cluster_hist_index = 2)
    ]

    area_spec = GhostedAreaSpec(ghost_maxrap = 1.0,
                                ghost_area = 0.05,
                                seed = 54321)

    csa = jet_reconstruct_area(particles;
                               algorithm = JetAlgorithm.AntiKt,
                               R = 0.4,
                               area_spec = area_spec,
                               strategy = RecoStrategy.N2Plain,)

    cs = cluster_sequence(csa)
    data = csa.area_data

    for i in (cs.n_initial_jets + 1):length(cs.history)
        history_element = cs.history[i]
        parent1 = history_element.parent1
        parent2 = history_element.parent2

        if parent2 == JetReconstruction.BeamJet
            @test data.areas[i] ≈ data.areas[parent1]
        elseif parent2 > 0
            @test data.areas[i] ≈ data.areas[parent1] + data.areas[parent2]
        end
    end
end

@testset "area calculation rejects ee algorithms" begin
    particles = [
        PseudoJet(100.0, 0.0, 0.0, 100.0; cluster_hist_index = 1)
    ]

    @test_throws ArgumentError jet_reconstruct_area(particles;
                                                    algorithm = JetAlgorithm.EEKt,
                                                    R = 0.4,)
end

@testset "area calculation rejects tiled strategies" begin
    particles = [
        PseudoJet(100.0, 0.0, 0.0, 100.0; cluster_hist_index = 1)
    ]

    @test_throws ArgumentError jet_reconstruct_area(particles;
                                                    algorithm = JetAlgorithm.AntiKt,
                                                    R = 0.4,
                                                    strategy = RecoStrategy.N2Tiled,)

    @test_throws ArgumentError jet_reconstruct_area(particles;
                                                    algorithm = JetAlgorithm.AntiKt,
                                                    R = 0.4,
                                                    strategy = RecoStrategy.Best,)
end