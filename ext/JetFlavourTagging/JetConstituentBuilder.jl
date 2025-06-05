module JetConstituentBuilder

using EDM4hep
using StructArrays: StructVector

# Define type aliases for clarity # TODO: Move to JetFlavourTagging file. But requires more checks.
const JetConstituents = StructVector{ReconstructedParticle}
const JetConstituentsData = Vector{Float32}

"""
    build_constituents(jets::JetConstituents, 
                       rps::Vector{EDM4hep.ReconstructedParticle}) -> Vector{JetConstituents}

Build the collection of constituents (mapping jet -> reconstructed particles) for all jets in event.

# Returns
A vector of JetConstituents, each containing the constituents for a specific jet.
"""
function build_constituents(jets::JetConstituents, 
                            rps::Vector{EDM4hep.ReconstructedParticle})
    jcs = Vector{JetConstituents}()
    for jet in jets
        constituents = JetConstituents()
        for i in jet.particles_begin:jet.particles_end-1
            push!(constituents, rps[i+1])  # Julia uses 1-based indexing
        end
        push!(jcs, constituents)
    end
    return jcs
end

"""
    build_constituents_cluster(rps::JetConstituents, 
                              indices::Vector{Vector{Int}}) -> Vector{JetConstituents}

Build the collection of constituents using cluster indices.

# Returns
A vector of JetConstituents, each containing the constituents for a specific cluster.
"""
function build_constituents_cluster(rps::JetConstituents, indices::Vector{Vector{Int}})
    jcs = Vector{JetConstituents}()
    for jet_indices in indices
        jc = rps[jet_indices]
        push!(jcs, jc)
    end
    return jcs
end
"""
    get_jet_constituents(csts::Vector{JetConstituents}, jet::Int) -> JetConstituents

Retrieve the constituents of an indexed jet in the event.

# Returns
The constituents of the specified jet, or an empty collection if the jet index is invalid.
"""
function get_jet_constituents(csts::Vector{JetConstituents}, jet::Int)
    if jet < 0
        return JetConstituents()
    end
    return csts[jet+1]  # Julia uses 1-based indexing
end

"""
    get_constituents(csts::Vector{JetConstituents}, jets::Vector{Int}) -> Vector{JetConstituents}

Retrieve the constituents of a collection of indexed jets in the event.
"""
function get_constituents(csts::Vector{JetConstituents}, jets::Vector{Int})
    jcs = Vector{JetConstituents}()
    for i in eachindex(jets)
        if jets[i] >= 0
            push!(jcs, csts[jets[i]+1])  # Julia uses 1-based indexing
        end
    end
    return jcs
end

end # module JetConstituentBuilder