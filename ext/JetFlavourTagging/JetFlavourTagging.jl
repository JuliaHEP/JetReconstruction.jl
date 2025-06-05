module JetFlavourTagging

using JetReconstruction
using EDM4hep
using JSON # TODO: Require to change this in main module
using ONNXRunTime
using StructArrays: StructVector
using LorentzVectorHEP


const JetConstituents = StructVector{ReconstructedParticle}
const JetConstituentsData = Vector{Float32}

# Building jet constituents from the EDM4hep event
include("JetConstituentBuilder.jl")

function JetReconstruction.build_constituents(jets::JetConstituents, 
                                                rps::Vector{EDM4hep.ReconstructedParticle})
    JetConstituentBuilder.build_constituents(jets, rps)
end

function JetReconstruction.build_constituents_cluster(rps::JetConstituents, 
                                                        indices::Vector{Vector{Int}})
    JetConstituentBuilder.build_constituents_cluster(rps, indices)
end

function JetReconstruction.get_jet_constituents(csts::Vector{JetConstituents}, jet::Int)
    JetConstituentBuilder.get_jet_constituents(csts, jet)
end

function JetReconstruction.get_constituents(csts::Vector{JetConstituents}, jets::Vector{Int})
    JetConstituentBuilder.get_constituents(csts, jets)
end

# Helper functions that talk with the ONNX model
include("JetFlavourHelper.jl")

# extract_features, setup_weaver, prepare_input_tensor, get_weights, get_weight
function JetReconstruction.extract_features(jets::Vector{EEJet}, jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                                            tracks::StructVector{EDM4hep.TrackState}, bz::Float32 = 2.0f0, 
                                            track_L::AbstractArray{T} where T <: AbstractFloat, 
                                            trackdata::StructVector{EDM4hep.Track}=StructVector{EDM4hep.Track}(), 
                                            trackerhits::StructVector{EDM4hep.TrackerHit}=StructVector{EDM4hep.TrackerHit}(), 
                                            gammadata::StructVector{EDM4hep.Cluster}=StructVector{EDM4hep.Cluster}(), 
                                            nhdata::StructVector{EDM4hep.Cluster}=StructVector{EDM4hep.Cluster}(), 
                                            calohits::StructVector{EDM4hep.CalorimeterHit}=StructVector{EDM4hep.CalorimeterHit}(), 
                                            dNdx::StructVector{EDM4hep.Quantity}=StructVector{EDM4hep.Quantity}())
    return JetFlavourHelper.extract_features(jets, jcs, tracks, bz, track_L, trackdata, trackerhits, gammadata, nhdata, calohits, dNdx)
end

function JetReconstruction.setup_weaver(onnx_path::String, json_path::String)
    return JetFlavourHelper.setup_weaver(onnx_path, json_path) # TODO: Change this to setup_onnx_runtime
end 

function JetReconstruction.prepare_input_tensor(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                            jets::Vector{EEJet}, 
                            config::Dict, 
                            feature_data::Dict)
    return JetFlavourHelper.prepare_input_tensor(jcs, jets, config, feature_data)
end

function JetReconstruction.get_weights(slot::Int, vars::Dict{String, Dict{String, Vector{Vector{Float32}}}}, 
                                        jets::Vector{EEJet}, jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                                        json_config::Dict, model::ONNXRunTime.InferenceSession)
    return JetFlavourHelper.get_weights(slot, vars, jets, jcs, json_config, model)
end

function JetReconstruction.get_weight(jet_weights::Vector{Vector{Float32}}, weight_idx::Int)
    return JetFlavourHelper.get_weight(jet_weights, weight_idx) 
end

end # module JetFlavourTagging