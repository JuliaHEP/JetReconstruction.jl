# Some basic tests of common utilities for jet types

include("common.jl")

pj1 = PseudoJet(1.0, 1.3, -2.3, 25.0)
pj2 = PseudoJet(1.0, 0.3, -1.3, 17.0)

eej1 = EEJet(1.0, 1.3, -2.3, 25.0)
eej2 = EEJet(-1.0, 3.2, -1.2, 39.0)

@testset "Common jet utilities" begin
    @test JetReconstruction.pt_fraction(pj1, pj2) ≈ 0.3889609897118418
    @test JetReconstruction.pt_fraction(eej1, eej2) ≈ 0.32850184248668773

    jr_kt_scale = JetReconstruction.kt_scale(pj1, pj2)
    lvhep_kt_scale = min(JetReconstruction.pt(pj1), JetReconstruction.pt(pj2)) *
                     deltar(JetReconstruction.lorentzvector_cyl(pj1),
                            JetReconstruction.lorentzvector_cyl(pj2))
    @test jr_kt_scale ≈ lvhep_kt_scale

    # Test conversions
    lv_pj1 = JetReconstruction.lorentzvector(pj1)
    lv_eej1 = JetReconstruction.lorentzvector(eej1)
    @test LorentzVectorHEP.energy(lv_pj1) ≈ JetReconstruction.energy(pj1) ≈
          JetReconstruction.energy(eej1) ≈ LorentzVectorHEP.energy(lv_eej1)
    @test LorentzVectorHEP.px(lv_pj1) ≈ JetReconstruction.px(pj1) ≈
          JetReconstruction.px(eej1) ≈ LorentzVectorHEP.px(lv_eej1)
    @test LorentzVectorHEP.py(lv_pj1) ≈ JetReconstruction.py(pj1) ≈
          JetReconstruction.py(eej1) ≈ LorentzVectorHEP.py(lv_eej1)
    @test LorentzVectorHEP.pz(lv_pj1) ≈ JetReconstruction.pz(pj1) ≈
          JetReconstruction.pz(eej1) ≈ LorentzVectorHEP.pz(lv_eej1)

    lvc_pj1 = JetReconstruction.lorentzvector_cyl(pj1)
    lvc_eej1 = JetReconstruction.lorentzvector_cyl(eej1)
    @test LorentzVectorHEP.pt(lvc_pj1) ≈ JetReconstruction.pt(eej1) ≈
          JetReconstruction.pt(eej1) ≈ LorentzVectorHEP.pt(lvc_eej1)
    @test LorentzVectorHEP.eta(lvc_pj1) ≈ JetReconstruction.eta(eej1) ≈
          JetReconstruction.eta(eej1) ≈ LorentzVectorHEP.eta(lvc_eej1)
    @test LorentzVectorHEP.phi(lvc_pj1) ≈ JetReconstruction.phi(eej1) ≈
          JetReconstruction.phi(eej1) ≈ LorentzVectorHEP.phi(lvc_eej1)
    @test LorentzVectorHEP.mass(lvc_pj1) ≈ JetReconstruction.mass(eej1) ≈
          JetReconstruction.mass(eej1) ≈ LorentzVectorHEP.mass(lvc_eej1)
end

@testset "Final jets extraction" begin
    # Test final jets extraction from different jet types
    vec_pj = [pj1, pj2]
    finaljets_pj = JetReconstruction.final_jets(vec_pj)

    vec_lorentzhep = [JetReconstruction.lorentzvector(pj1),
        JetReconstruction.lorentzvector(pj2)]
    finaljets_lorentzhep = JetReconstruction.final_jets(vec_lorentzhep)

    vec_lorentzhepcyl = [JetReconstruction.lorentzvector_cyl(pj1),
        JetReconstruction.lorentzvector_cyl(pj2)]
    finaljets_lorentzhepcyl = JetReconstruction.final_jets(vec_lorentzhepcyl)

    @test length(finaljets_pj) == 2
    @test length(finaljets_lorentzhep) == 2
    @test length(finaljets_lorentzhepcyl) == 2
    for i in 1:2
        @test finaljets_pj[i] ≈ finaljets_lorentzhep[i] ≈ finaljets_lorentzhepcyl[i]
    end
end
