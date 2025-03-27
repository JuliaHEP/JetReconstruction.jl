# Some basic tests of common utilities for jet types

include("common.jl")

pj1 = PseudoJet(1.0, 1.3, -2.3, 25.0)
pj2 = PseudoJet(1.0, 0.3, -1.3, 17.0)

eej1 = EEjet(1.0, 1.2, -2.2, 22.0)
eej2 = EEjet(-1.0, 3.2, -1.2, 39.0)

@testset "Common jet utilities" begin
    @test JetReconstruction.pt_fraction(pj1, pj2) ≈ 0.3889609897118418
    @test JetReconstruction.pt_fraction(eej1, eej2) ≈ 0.3178347357639779

    jr_kt_scale = JetReconstruction.kt_scale(pj1, pj2)
    lvhep_kt_scale = min(JetReconstruction.pt(pj1), JetReconstruction.pt(pj2)) *
                     deltar(JetReconstruction.lorentzvector_cyl(pj1),
                            JetReconstruction.lorentzvector_cyl(pj2))
    @test jr_kt_scale ≈ lvhep_kt_scale

    # Test conversions
    lv1 = JetReconstruction.lorentzvector(pj1)
    @test LorentzVectorHEP.energy(lv1) ≈ JetReconstruction.energy(pj1)
    @test LorentzVectorHEP.px(lv1) ≈ JetReconstruction.px(pj1)
    @test LorentzVectorHEP.py(lv1) ≈ JetReconstruction.py(pj1)
    @test LorentzVectorHEP.pz(lv1) ≈ JetReconstruction.pz(pj1)

    lv_cyl1 = JetReconstruction.lorentzvector_cyl(eej1)
    @test LorentzVectorHEP.pt(lv_cyl1) ≈ JetReconstruction.pt(eej1)
    @test LorentzVectorHEP.eta(lv_cyl1) ≈ JetReconstruction.eta(eej1)
    @test LorentzVectorHEP.phi(lv_cyl1) ≈ JetReconstruction.phi(eej1)
    @test LorentzVectorHEP.mass(lv_cyl1) ≈ JetReconstruction.mass(eej1)
end
