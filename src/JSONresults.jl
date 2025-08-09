"""
    struct FinalJet

A struct representing the final properties of a jet, used for
JSON serialisation.

# Fields
- `rap::Float64`: The rapidity of the jet.
- `phi::Float64`: The azimuthal angle of the jet.
- `pt::Float64`: The transverse momentum of the jet.
"""
struct FinalJet
    rap::Float64
    phi::Float64
    pt::Float64
end

"""
    FinalJet(jet::T)

Convert a jet of type `T` to a `FinalJet` object. `T` must be a type that
supports the methods `rapidity`, `phi`, and `pt`.
"""
FinalJet(jet::T) where {T} = FinalJet(rapidity(jet), phi(jet), pt(jet))

import Base: ==, ≈
"""
    ==(fj1::FinalJet, fj2::FinalJet)

Compare two `FinalJet` objects for equality based on their rapidity, azimuthal
angle, and transverse momentum.
"""
function ==(fj1::FinalJet, fj2::FinalJet)
    (fj1.rap == fj2.rap) && (fj1.phi == fj2.phi) && (fj1.pt == fj2.pt)
end

"""
    ≈(fj1::FinalJet, fj2::FinalJet)

Compare two `FinalJet` objects for approximate equality based on their rapidity,
azimuthal angle, and transverse momentum.
"""
function ≈(fj1::FinalJet, fj2::FinalJet; tol = 1e-6)
    fj1.rap ≈ fj2.rap && fj1.phi ≈ fj2.phi && fj1.pt ≈ fj2.pt
end

"""
    struct FinalJets

A struct with the vector of all jets for a certain jet identifier, used for
JSON serialisation.

# Fields
- `jetid::Int64`: The ID of the jet.
- `jets::Vector{FinalJet}`: A vector of `FinalJet` objects representing the jets.
"""
struct FinalJets
    jetid::Int64
    jets::Vector{FinalJet}
end
