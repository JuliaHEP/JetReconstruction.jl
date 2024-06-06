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
