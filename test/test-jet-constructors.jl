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
    @test_throws ArgumentError _=PseudoJet{Float32}(broken4v)
    @test_throws ArgumentError _=EEJet{Float32}(broken4v)

    pj_generic_typed = PseudoJet{Float32}(g4v)
    @test typeof(JetReconstruction.px(pj_generic_typed)) == Float32
    eej_generic_typed = EEJet{Float32}(g4v)
    @test typeof(JetReconstruction.px(eej_generic_typed)) == Float32
end

# Test construction from LorentzVector and LorentzVectorCyl
lv = LorentzVector(25.0, 1.0, 1.3, -2.3)
lv_cyl = LorentzVectorCyl(pt(lv), eta(lv), phi(lv), mass(lv))
@testset "LorentzVector construction" begin
    pj_from_lv = PseudoJet(lv)
    @test JetReconstruction.px(pj_from_lv) == px(lv)
    @test JetReconstruction.py(pj_from_lv) == py(lv)
    @test JetReconstruction.pz(pj_from_lv) == pz(lv)
    @test JetReconstruction.energy(pj_from_lv) == energy(lv)

    eej_from_lv = EEJet(lv)
    @test JetReconstruction.px(eej_from_lv) == px(lv)
    @test JetReconstruction.py(eej_from_lv) == py(lv)
    @test JetReconstruction.pz(eej_from_lv) == pz(lv)
    @test JetReconstruction.energy(eej_from_lv) == energy(lv)

    eej_typed_from_lv = EEJet{Float32}(lv)
    @test typeof(JetReconstruction.px(eej_typed_from_lv)) == Float32
    pj_typed_from_lv = PseudoJet{Float32}(lv)
    @test typeof(JetReconstruction.px(pj_typed_from_lv)) == Float32
end

@testset "LorentzVectorCyl construction" begin
    pj_from_lvc = PseudoJet(lv_cyl)
    @test JetReconstruction.px(pj_from_lvc) ≈ px(lv_cyl)
    @test JetReconstruction.py(pj_from_lvc) ≈ py(lv_cyl)
    @test JetReconstruction.pz(pj_from_lvc) ≈ pz(lv_cyl)
    @test JetReconstruction.energy(pj_from_lvc) ≈ energy(lv_cyl)

    eej_from_lvc = EEJet(lv_cyl)
    @test JetReconstruction.px(eej_from_lvc) ≈ px(lv_cyl)
    @test JetReconstruction.py(eej_from_lvc) ≈ py(lv_cyl)
    @test JetReconstruction.pz(eej_from_lvc) ≈ pz(lv_cyl)
    @test JetReconstruction.energy(eej_from_lvc) ≈ energy(lv_cyl)

    eej_typed_from_lvc = EEJet{Float32}(lv_cyl)
    @test typeof(JetReconstruction.px(eej_typed_from_lvc)) == Float32
    pj_typed_from_lvc = PseudoJet{Float32}(lv_cyl)
    @test typeof(JetReconstruction.px(pj_typed_from_lvc)) == Float32
end

# Mixed type constructors - check that promotion works as expected
@testset "Mixed type constructors and conversions" begin
    eej_promote = EEJet(1.0f0, 2.0f0, 3.0f0, 40.0)
    @test typeof(JetReconstruction.px(eej_promote)) == Float64
    pj_promote = PseudoJet(1.0f0, 2.0f0, 3.0f0, 40.0)
    @test typeof(JetReconstruction.px(pj_promote)) == Float64

    eej_demote = EEJet{Float32}(eej_promote)
    @test typeof(JetReconstruction.px(eej_demote)) == Float32
    pj_demote = PseudoJet{Float32}(pj_promote)
    @test typeof(JetReconstruction.px(pj_demote)) == Float32
end
