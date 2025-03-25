# Some basic tests of common utilities for jet types

include("common.jl")

toLVHEP(jet) = fromPtEtaPhiE(JetReconstruction.pt(jet), JetReconstruction.eta(jet), JetReconstruction.phi(jet), JetReconstruction.energy(jet))

pj1=PseudoJet(1.0,1.3,-2.3,25.0)
pj2=PseudoJet(1.0,0.3,-1.3,17.0)

# lvhep1=toLVHEP(pj1)
# lvhep2=toLVHEP(pj2)

eej1=EEjet(1.0,1.2,-2.2,22.0)
eej2=EEjet(-1.0,3.2,-1.2,39.0)

@testset "Common jet utilities" begin
    @test JetReconstruction.momentum_fraction(pj1,pj2) ≈ 0.3889609897118418
    @test JetReconstruction.momentum_fraction(eej1,eej2) ≈ 0.3178347357639779

    jr_kt_scale = JetReconstruction.kt_scale(pj1,pj2)
    lvhep_kt_scale = min(JetReconstruction.pt(pj1), JetReconstruction.pt(pj2)) * deltar(toLVHEP(pj1), toLVHEP(pj2))
    @test jr_kt_scale ≈ lvhep_kt_scale
end
