module JetFlavourTagging

using JetReconstruction
using EDM4hep
using JSON
using ONNXRunTime
using StructArrays: StructVector
using LorentzVectorHEP

const JetConstituents = StructVector{ReconstructedParticle, <:Any}
const JetConstituentsData = Vector{Float32}

# Building jet constituents from the EDM4hep event
include("JetConstituentBuilder.jl")

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
function JetReconstruction.build_constituents_cluster(reco_particles::JetConstituents,
                                                      indices::Vector{Vector{Int64}})
    JetConstituentBuilder.build_constituents_cluster(reco_particles, indices)
end

# TODO: Restore get_jets_constituents, get_constituents, build_constituents function.

# Helper functions that talk with the ONNX model
include("JetFlavourHelper.jl")

"""
    extract_features(jets::Vector{EEJet}, jets_constituents::Vector{<:JetConstituents}, 
                    tracks::AbstractVector{EDM4hep.TrackState}, bz::Float32, 
                    track_L::AbstractArray{T} where T <: AbstractFloat,
                    json_config::Dict,
                    trackdata::AbstractVector{EDM4hep.Track}=AbstractVector{EDM4hep.Track}(), 
                    trackerhits::AbstractVector{EDM4hep.TrackerHit}=AbstractVector{EDM4hep.TrackerHit}(), 
                    gammadata::AbstractVector{EDM4hep.Cluster}=AbstractVector{EDM4hep.Cluster}(), 
                    nhdata::AbstractVector{EDM4hep.Cluster}=AbstractVector{EDM4hep.Cluster}(), 
                    calohits::AbstractVector{EDM4hep.CalorimeterHit}=AbstractVector{EDM4hep.CalorimeterHit}(), 
                    dNdx::AbstractVector{EDM4hep.Quantity}=AbstractVector{EDM4hep.Quantity}()) -> Dict

Extract features for jet flavour tagging based on JSON configuration.

# Arguments 
- `jets`: Vector of jets (EEJet)
- `jets_constituents`: Vector of jet constituents
- `tracks`: StructVector of track states
- `bz`: Magnetic field strength
- `track_L`: Array of track lengths
- `json_config`: JSON configuration dict specifying which features to extract (required)
- `trackdata`: Vector of track data (optional)
- `trackerhits`: Vector of tracker hits (optional)
- `gammadata`: Vector of gamma clusters (optional)
- `nhdata`: Vector of neutral hadron clusters (optional)
- `calohits`: Vector of calorimeter hits (optional)
- `dNdx`: Vector of dE/dx measurements (optional)   

# Returns
Dictionary containing extracted features as specified in the JSON configuration.
"""
function JetReconstruction.extract_features(jets::Vector{EEJet},
                                            jets_constituents::Vector{<:JetConstituents},
                                            tracks::AbstractVector{EDM4hep.TrackState},
                                            bz::Float32,
                                            track_L::AbstractArray{T} where {T <:
                                                                             AbstractFloat},
                                            json_config::Dict,
                                            trackdata::AbstractVector{EDM4hep.Track} = AbstractVector{EDM4hep.Track}(),
                                            trackerhits::AbstractVector{EDM4hep.TrackerHit} = AbstractVector{EDM4hep.TrackerHit}(),
                                            gammadata::AbstractVector{EDM4hep.Cluster} = AbstractVector{EDM4hep.Cluster}(),
                                            nhdata::AbstractVector{EDM4hep.Cluster} = AbstractVector{EDM4hep.Cluster}(),
                                            calohits::AbstractVector{EDM4hep.CalorimeterHit} = AbstractVector{EDM4hep.CalorimeterHit}(),
                                            dNdx::AbstractVector{EDM4hep.Quantity} = AbstractVector{EDM4hep.Quantity}())
    return JetFlavourHelper.extract_features(jets, jets_constituents, tracks, bz, track_L,
                                             json_config, trackdata,
                                             trackerhits, gammadata, nhdata, calohits, dNdx)
end

"""
    setup_onnx_runtime(onnx_path::AbstractString, json_path::AbstractString) -> ONNXRunTime.InferenceSession

Setup the ONNX model and preprocessing configuration for jet flavour tagging.

# Arguments
- `onnx_path`: Path to the ONNX model file
- `json_path`: Path to the JSON configuration file

# Returns
An ONNX inference session for the loaded model
"""
function JetReconstruction.setup_onnx_runtime(onnx_path::AbstractString,
                                              json_path::AbstractString)
    return JetFlavourHelper.setup_onnx_runtime(onnx_path, json_path)
end

"""
    prepare_input_tensor(jets_constituents::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                        jets::Vector{EEJet}, 
                        config::Dict, 
                        feature_data::Dict) -> Dict{String, Array}

Prepare input tensors for the neural network from jet constituents.

# Arguments
- `jets_constituents`: Vector of jet constituents (structured as a vector of StructVector of ReconstructedParticle)
- `jets`: Vector of jets (EEJet)
- `config`: JSON configuration for preprocessing
- `feature_data`: Dictionary containing all extracted features

# Returns
Dictionary of input tensors
"""
function JetReconstruction.prepare_input_tensor(jets_constituents::Vector{<:JetConstituents},
                                                jets::Vector{EEJet},
                                                config::Dict,
                                                feature_data::Dict)
    return JetFlavourHelper.prepare_input_tensor(jets_constituents, jets, config,
                                                 feature_data)
end

"""
    get_weights(slot::Int, vars::Dict{String, Vector{Vector{Float32}}}, 
                jets::Vector{EEJet}, json_config::Dict, model::ONNXRunTime.InferenceSession) -> Vector{Vector{Float32}}

Compute jet flavour probabilities for each jet.

# Arguments
- `slot`: Threading slot
- `vars`: Dictionary containing all features for jet constituents
- `jets`: Vector of jets
- `json_config`: JSON configuration for preprocessing
- `model`: ONNX inference session

# Returns
Vector of flavor probabilities for each jet
"""
function JetReconstruction.get_weights(slot::Int,
                                       vars::Dict{String,
                                                  Dict{String, Vector{Vector{Float32}}}},
                                       jets::Vector{EEJet},
                                       jets_constituents::Vector{<:JetConstituents},
                                       json_config::Dict,
                                       model::ONNXRunTime.InferenceSession)
    return JetFlavourHelper.get_weights(slot, vars, jets, jets_constituents, json_config,
                                        model)
end

"""
    get_weight(jet_weights::Vector{Vector{Float32}}, weight_idx::Int) -> Vector{Float32}

Extract a specific weight/score from the jet weights.

# Arguments
- `jet_weights`: Vector of weight vectors for each jet
- `weight_idx`: Index of the weight to extract

# Returns
Vector of the specified weight for each jet
"""
function JetReconstruction.get_weight(jet_weights::Vector{Vector{Float32}}, weight_idx::Int)
    return JetFlavourHelper.get_weight(jet_weights, weight_idx)
end

end # module JetFlavourTagging
