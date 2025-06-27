module JetConstituentBuilder

using EDM4hep
using StructArrays: StructVector

# Define type aliases for clarity # TODO: Move to JetFlavourTagging file. But requires more checks.
const JetConstituents = StructVector{ReconstructedParticle, <:Any}
const JetConstituentsData = Vector{Float32}

"""
    build_constituents(jets::JetConstituents, 
                       rps::JetConstituents) -> Vector{JetConstituents}

Build the collection of constituents (mapping jet -> reconstructed particles) for all jets in event.

# Returns
A vector of JetConstituents, each containing the constituents for a specific jet.
"""
# TODO: Fix this function to be interpolate with Julia pipeline. Specificly, what would be the input jets?

"""
    build_constituents_cluster(rps::JetConstituents, 
                               indices::Vector{Vector{Int}}) -> Vector{JetConstituents}

Build the collection of constituents using cluster indices.

# Arguments
- rps: a vector of `JetConstituents` representing reconstructed particles.
- indices: a vector of vectors, where each inner vector contains indices of particles for a specific cluster.

# Returns
A vector of JetConstituents, each containing the constituents for a specific cluster.
"""
function build_constituents_cluster(rps, indices)
    return map(jet_indices -> rps[jet_indices], indices)
end

"""
    get_jet_constituents(csts::Vector{JetConstituents}, jet_index::Int) -> JetConstituents

Retrieve the constituents of an indexed jet in the event.
# Arguments
- csts: constituents collection, a vector of `JetConstituents`.
- jet_index: the index of the jet for which to retrieve constituents (1-based index) 

# Returns
The constituents of the specified jet, or an empty collection if the jet index is invalid.
"""
function get_jet_constituents(csts::Vector{JetConstituents}, jet_index::Int)
    return csts[jet_index]
end

"""
    get_constituents(csts::Vector{JetConstituents}, jet_indices::Vector{Int}) -> Vector{JetConstituents}

Retrieve the constituents of a collection of indexed jets in the event.

# Arguments
- csts: constituents collection, a vector of `JetConstituents`.
- jet_indices: a vector of jet indices (1-based index) for which to retrieve constituents.

# Returns
A vector of `JetConstituents`, each containing the constituents for the specified jets.
"""
function get_constituents(csts::Vector{JetConstituents}, jet_indices::Vector{Int})
    jcs = Vector{JetConstituents}()
    for i in 1:length(jet_indices)
        if jet_indices[i] >= 1
            if i <= length(csts)
                push!(jcs, csts[i])
            end
        end
    end
    return jcs
end

end # module JetConstituentBuilder
