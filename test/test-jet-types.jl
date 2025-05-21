# A few very basic tests for the Jet types.

include("common.jl")

pj = PseudoJet(1.0, 2.0, 3.0, 10.0; cluster_hist_index = 7)
pj_beam = PseudoJet(0.0, 0.0, 5.0, 5.0)
@testset "PseudoJet tests" begin
    @test isbits(pj)
    @test JetReconstruction.px(pj) ≈ 1.0
    @test JetReconstruction.py(pj) ≈ 2.0
    @test JetReconstruction.pz(pj) ≈ 3.0
    @test JetReconstruction.energy(pj) ≈ 10.0
    @test JetReconstruction.mag(pj) ≈ sqrt(1.0^2 + 2.0^2 + 3.0^2)
    @test JetReconstruction.pt2(pj) ≈ 1.0^2 + 2.0^2
    @test JetReconstruction.phi(pj) ≈ atan(2.0, 1.0)
    @test JetReconstruction.mass(pj) ≈ sqrt(10.0^2 - 1.0^2 - 2.0^2 - 3.0^2)
    @test JetReconstruction.cluster_hist_index(pj) == 7
    @test isvalid(pj) == true

    @test JetReconstruction.rapidity(pj_beam) ≈
          JetReconstruction._MaxRap + JetReconstruction.pz(pj_beam)
    @test JetReconstruction.eta(pj_beam) ≈ JetReconstruction._MaxRap

    # This isn't really a test of the output, but rather that the object
    # can be printed without error
    @test string(pj) == "PseudoJet(px: 1.0 py: 2.0 pz: 3.0 E: 10.0 cluster_hist_index: 7)"
end

eej = EEJet(1.0, 2.0, 3.0, 10.0; cluster_hist_index = 7)
eej_beam = EEJet(0.0, 0.0, 5.0, 5.0)
@testset "EEJet tests" begin
    @test isbits(eej)
    @test JetReconstruction.px(eej) ≈ 1.0
    @test JetReconstruction.py(eej) ≈ 2.0
    @test JetReconstruction.pz(eej) ≈ 3.0
    @test JetReconstruction.energy(eej) ≈ 10.0
    @test JetReconstruction.p2(eej) ≈ 1.0^2 + 2.0^2 + 3.0^2
    @test JetReconstruction.pt2(eej) ≈ 1.0^2 + 2.0^2
    @test JetReconstruction.phi(eej) ≈ atan(2.0, 1.0)
    @test JetReconstruction.mass(eej) ≈ sqrt(10.0^2 - 1.0^2 - 2.0^2 - 3.0^2)
    @test JetReconstruction.cluster_hist_index(eej) == 7

    @test JetReconstruction.rapidity(eej_beam) ≈
          JetReconstruction._MaxRap + JetReconstruction.pz(eej_beam)

    # This isn't really a test of the output, but rather that the object
    # can be printed without error
    @test string(eej) == "EEJet(px: 1.0 py: 2.0 pz: 3.0 E: 10.0 cluster_hist_index: 7)"
end

eej_pseudojet = EEJet(pj)
pseudojet_eej = PseudoJet(eej)
@testset "Jet Interoperability" begin
    @test JetReconstruction.px(eej_pseudojet) ≈ JetReconstruction.px(pj)
    @test JetReconstruction.py(eej_pseudojet) ≈ JetReconstruction.py(pj)
    @test JetReconstruction.pz(eej_pseudojet) ≈ JetReconstruction.pz(pj)
    @test JetReconstruction.energy(eej_pseudojet) ≈ JetReconstruction.energy(pj)

    @test JetReconstruction.px(pseudojet_eej) ≈ JetReconstruction.px(eej)
    @test JetReconstruction.py(pseudojet_eej) ≈ JetReconstruction.py(eej)
    @test JetReconstruction.pz(pseudojet_eej) ≈ JetReconstruction.pz(eej)
    @test JetReconstruction.energy(pseudojet_eej) ≈ JetReconstruction.energy(eej)
end
