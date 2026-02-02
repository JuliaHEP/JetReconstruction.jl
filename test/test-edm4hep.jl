# Test that we can work with EDM4hep particles

include("common.jl")

using EDM4hep

dummyRecoParticle = ReconstructedParticle(energy = 4.0f0,
                                          momentum = Vector3f(1.0f0, 2.0f0, 3.0f0))

@testset "Construction of jets from EDM4hep particles" begin
    eej = EEJet(dummyRecoParticle)
    pj = PseudoJet(dummyRecoParticle)
    @test typeof(EEJet(dummyRecoParticle)) === EEJet
    @test typeof(EEJet(dummyRecoParticle; cluster_hist_index = 99)) === EEJet
    @test JetReconstruction.energy(eej) == dummyRecoParticle.energy
    @test JetReconstruction.px(eej) == dummyRecoParticle.momentum.x
    @test JetReconstruction.py(eej) == dummyRecoParticle.momentum.y
    @test JetReconstruction.pz(eej) == dummyRecoParticle.momentum.z

    @test typeof(PseudoJet(dummyRecoParticle)) === PseudoJet
    @test typeof(PseudoJet(dummyRecoParticle; cluster_hist_index = 99)) === PseudoJet
    @test JetReconstruction.energy(pj) == dummyRecoParticle.energy
    @test JetReconstruction.px(pj) == dummyRecoParticle.momentum.x
    @test JetReconstruction.py(pj) == dummyRecoParticle.momentum.y
    @test JetReconstruction.pz(pj) == dummyRecoParticle.momentum.z
end
