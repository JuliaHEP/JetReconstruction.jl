# A few very basic tests for the Jet types.

include("common.jl")

pj = PseudoJet(1.0, 2.0, 3.0, 10.0)
pj._cluster_hist_index = 7
@testset "PseudoJet tests" begin
    @test JetReconstruction.px(pj) ≈ 1.0
    @test JetReconstruction.py(pj) ≈ 2.0
    @test JetReconstruction.pz(pj) ≈ 3.0
    @test JetReconstruction.energy(pj) ≈ 10.0
    @test JetReconstruction.mag(pj) ≈ sqrt(1.0^2 + 2.0^2 + 3.0^2)
    @test JetReconstruction.pt2(pj) ≈ 1.0^2 + 2.0^2
    @test JetReconstruction.phi(pj) ≈ atan(2.0, 1.0)
    @test JetReconstruction.mass(pj) ≈ sqrt(10.0^2 - 1.0^2 - 2.0^2 - 3.0^2)
    @test JetReconstruction.cluster_hist_index(pj) == 7
end

eej = EEJet(1.0, 2.0, 3.0, 10.0, 7)
@testset "EEJet tests" begin
    @test JetReconstruction.px(eej) ≈ 1.0
    @test JetReconstruction.py(eej) ≈ 2.0
    @test JetReconstruction.pz(eej) ≈ 3.0
    @test JetReconstruction.energy(eej) ≈ 10.0
    @test JetReconstruction.p2(eej) ≈ 1.0^2 + 2.0^2 + 3.0^2
    @test JetReconstruction.pt2(eej) ≈ 1.0^2 + 2.0^2
    @test JetReconstruction.phi(eej) ≈ atan(2.0, 1.0)
    @test JetReconstruction.mass(eej) ≈ sqrt(10.0^2 - 1.0^2 - 2.0^2 - 3.0^2)
    @test JetReconstruction.cluster_hist_index(eej) == 7
end

eej_pseudojet = EEJet(pj)
@testset "Jet Interoperability" begin
    @test JetReconstruction.px(eej_pseudojet) ≈ JetReconstruction.px(pj)
    @test JetReconstruction.py(eej_pseudojet) ≈ JetReconstruction.py(pj)
    @test JetReconstruction.pz(eej_pseudojet) ≈ JetReconstruction.pz(pj)
    @test JetReconstruction.energy(eej_pseudojet) ≈ JetReconstruction.energy(pj)
end
