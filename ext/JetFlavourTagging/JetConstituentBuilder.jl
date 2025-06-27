module JetConstituentBuilder

using EDM4hep
using StructArrays: StructVector

# Define type aliases for clarity
const JetConstituents = StructVector{ReconstructedParticle, <:Any}
const JetConstituentsData = Vector{Float32}

"""
    build_constituents(jets::JetConstituents, 
                       reco_particles::JetConstituents) -> Vector{JetConstituents}

Build the collection of constituents (mapping jet -> reconstructed particles) for all jets in event.

# Returns
A vector of JetConstituents, each containing the constituents for a specific jet.
"""
# TODO: Fix this function to be interpolate with Julia pipeline. Specificly, what would be the input jets?

"""
    build_constituents_cluster(reco_particles::JetConstituents, 
                               indices::Vector{Vector{Int}}) -> Vector{JetConstituents}

Build the collection of constituents using cluster indices.

# Arguments
- reco_particles: a vector of `JetConstituents` representing reconstructed particles.
- indices: a vector of vectors, where each inner vector contains indices of particles for a specific cluster.

# Returns
A vector of JetConstituents, each containing the constituents for a specific cluster.
"""
function build_constituents_cluster(reco_particles, indices)
    return map(jet_indices -> reco_particles[jet_indices], indices)
end

"""
    get_jet_constituents(constituents_collection::Vector{JetConstituents}, jet_index::Int) -> JetConstituents

Retrieve the constituents of an indexed jet in the event.
# Arguments
- constituents_collection: constituents collection, a vector of `JetConstituents`.
- jet_index: the index of the jet for which to retrieve constituents (1-based index) 

# Returns
The constituents of the specified jet, or an empty collection if the jet index is invalid.
"""
function get_jet_constituents(constituents_collection::Vector{JetConstituents},
                              jet_index::Int)
    return constituents_collection[jet_index]
end

"""
    get_constituents(constituents_collection::Vector{JetConstituents}, jet_indices::Vector{Int}) -> Vector{JetConstituents}

Retrieve the constituents of a collection of indexed jets in the event.

# Arguments
- constituents_collection: constituents collection, a vector of `JetConstituents`.
- jet_indices: a vector of jet indices (1-based index) for which to retrieve constituents.

# Returns
A vector of `JetConstituents`, each containing the constituents for the specified jets.
"""
function get_constituents(constituents_collection::Vector{JetConstituents},
                          jet_indices::Vector{Int})
    # Filter valid indices and map to corresponding constituents
    valid_indices = filter(idx -> 1 <= idx <= length(constituents_collection), jet_indices)
    return map(idx -> constituents_collection[idx], valid_indices)
end

end # module JetConstituentBuilder
