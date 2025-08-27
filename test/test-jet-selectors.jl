# Test that final jets can be retrieved with different
# output types

include("common.jl")

"""
    test_jet_equality(ref, test)

Test two jet representations for equality.
"""
function test_jet_equality(ref, test)
    for (j1, j2) in zip(ref, test)
        @test JetReconstruction.pt(j1) ≈ JetReconstruction.pt(j2)
        test_phi = JetReconstruction.phi02pi(j2)
        test_phi = test_phi < 0.0 ? test_phi + 2π : test_phi
        @test JetReconstruction.phi(j1) ≈ test_phi
        @test JetReconstruction.rapidity(j1) ≈ JetReconstruction.rapidity(j2)
        @test JetReconstruction.eta(j1) ≈ JetReconstruction.eta(j2)
        @test JetReconstruction.px(j1) ≈ JetReconstruction.px(j2)
        @test JetReconstruction.py(j1) ≈ JetReconstruction.py(j2)
        @test JetReconstruction.pz(j1) ≈ JetReconstruction.pz(j2)
        @test JetReconstruction.energy(j1) ≈ JetReconstruction.energy(j2)
        @test JetReconstruction.mass2(j1) ≈ JetReconstruction.mass2(j2)
        @test JetReconstruction.mass(j1) ≈ JetReconstruction.mass(j2)
    end
    nothing
end

@testset "Jet selections to types" begin
    pp_event = first(JetReconstruction.read_final_state_particles(events_file_pp;
                                                                  maxevents = 1))
    ee_event = first(JetReconstruction.read_final_state_particles(events_file_ee;
                                                                  maxevents = 1))

    pp_cs = jet_reconstruct(pp_event; algorithm = JetAlgorithm.Kt)
    ee_cs = jet_reconstruct(ee_event; algorithm = JetAlgorithm.Durham)

    pp_ref_jets = inclusive_jets(pp_cs, PseudoJet; ptmin = 10.0)
    for T in (LorentzVectorCyl, LorentzVector)
        test_jet_equality(pp_ref_jets, inclusive_jets(pp_cs, T; ptmin = 10.0))
    end

    ee_ref_jets = exclusive_jets(ee_cs, EEJet; njets = 2)
    for T in (LorentzVectorCyl, LorentzVector)
        test_jet_equality(ee_ref_jets, exclusive_jets(ee_cs, T; njets = 2))
    end

    # Also check that for an unknown type we throw
    @test_throws ErrorException begin
        inclusive_jets(pp_cs, Float64; ptmin = 10.0)
    end
    @test_throws ErrorException begin
        exclusive_jets(ee_cs, Float64; njets = 2)
    end
end
