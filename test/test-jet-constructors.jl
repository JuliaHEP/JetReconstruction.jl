# Test that our internal jet types can be correctly constructed
# from the supported input 4-momentum types
include("common.jl")

# Define a different custom jet type that implements the LorentzVectorBase
# interface to test that our jet types can be constructed from these
struct genericPxPyPzE
    xmom::Real
    ymom::Real
    zmom::Real
    energy::Real
end
LorentzVectorBase.coordinate_system(::genericPxPyPzE) = LorentzVectorBase.PxPyPzE()
LorentzVectorBase.px(jet::genericPxPyPzE) = jet.xmom
LorentzVectorBase.py(jet::genericPxPyPzE) = jet.ymom
LorentzVectorBase.pz(jet::genericPxPyPzE) = jet.zmom
LorentzVectorBase.E(jet::genericPxPyPzE) = jet.energy

struct willnotwork
    xmom::Real
    ymom::Real
    zmom::Real
    energy::Real
end

g4v = genericPxPyPzE(1.0, 1.3, -2.3, 25.0)
broken4v = willnotwork(1.0, 1.3, -2.3, 25.0)

@testset "LorentzVectorBase interface construction" begin
    pj_generic = PseudoJet(g4v)
    @test JetReconstruction.px(pj_generic) == g4v.xmom
    @test JetReconstruction.py(pj_generic) == g4v.ymom
    @test JetReconstruction.pz(pj_generic) == g4v.zmom
    @test JetReconstruction.energy(pj_generic) == g4v.energy

    ej_generic = EEJet(g4v)
    @test JetReconstruction.px(ej_generic) == g4v.xmom
    @test JetReconstruction.py(ej_generic) == g4v.ymom
    @test JetReconstruction.pz(ej_generic) == g4v.zmom
    @test JetReconstruction.energy(ej_generic) == g4v.energy

    @test_throws ArgumentError _=PseudoJet(broken4v)
    @test_throws ArgumentError _=EEJet(broken4v)
end

# Test construction from LorentzVector and LorentzVectorCyl
lv = LorentzVector(25.0, 1.0, 1.3, -2.3)
lv_cyl = LorentzVectorCyl(pt(lv), eta(lv), phi(lv), mass(lv))
@testset "LorentzVector construction" begin
    eej_from_lv = PseudoJet(lv)
    @test JetReconstruction.px(eej_from_lv) == px(lv)
    @test JetReconstruction.py(eej_from_lv) == py(lv)
    @test JetReconstruction.pz(eej_from_lv) == pz(lv)
    @test JetReconstruction.energy(eej_from_lv) == energy(lv)

    eej_from_lv = EEJet(lv)
    @test JetReconstruction.px(eej_from_lv) == px(lv)
    @test JetReconstruction.py(eej_from_lv) == py(lv)
    @test JetReconstruction.pz(eej_from_lv) == pz(lv)
    @test JetReconstruction.energy(eej_from_lv) == energy(lv)
end

@testset "LorentzVectorCyl construction" begin
    eej_from_lv = PseudoJet(lv_cyl)
    @test JetReconstruction.px(eej_from_lv) ≈ px(lv_cyl)
    @test JetReconstruction.py(eej_from_lv) ≈ py(lv_cyl)
    @test JetReconstruction.pz(eej_from_lv) ≈ pz(lv_cyl)
    @test JetReconstruction.energy(eej_from_lv) ≈ energy(lv_cyl)

    eej_from_lv = EEJet(lv)
    @test JetReconstruction.px(eej_from_lv) ≈ px(lv_cyl)
    @test JetReconstruction.py(eej_from_lv) ≈ py(lv_cyl)
    @test JetReconstruction.pz(eej_from_lv) ≈ pz(lv_cyl)
    @test JetReconstruction.energy(eej_from_lv) ≈ energy(lv_cyl)
end
