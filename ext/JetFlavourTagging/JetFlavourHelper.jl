module JetFlavourHelper

using JSON
using ONNXRunTime
using StructArrays: StructVector
using JetReconstruction
using EDM4hep
using LorentzVectorHEP

# Utility functions for the jet features from EDM4hep
include("JetConstituentUtils.jl")
# include("ReconstructedParticle.jl")

"""
    JetFlavourHelper

A module for jet flavour identification using neural networks.
"""

"""
    setup_onnx_runtime(onnx_path::AbstractString, json_path::AbstractString) -> ONNXRunTime.InferenceSession

Setup the ONNX model and preprocessing configuration for jet flavour tagging.

# Arguments
- `onnx_path`: Path to the ONNX model file
- `json_path`: Path to the JSON configuration file

# Returns
An ONNX inference session for the loaded model
"""
function setup_onnx_runtime(onnx_path::AbstractString, json_path::AbstractString)
    # Load JSON configuration
    config = JSON.parsefile(json_path)
    model = ONNXRunTime.load_inference(onnx_path)

    return model, config
end

"""
    normalize_feature(value::Float32, info::Dict) -> Float32

Normalize a feature value based on the preprocessing information.

# Arguments
- `value`: Raw feature value
- `info`: Dictionary containing normalization parameters

# Returns
Normalized feature value
"""
function normalize_feature(value::Float32, info::Dict)
    if value == -9.0f0
        return 0.0f0  # Replace -9.0 (missing value) with 0
    end

    # Apply normalization using median and norm_factor
    normalized = (value - info["median"]) * info["norm_factor"]

    # Clamp to specified bounds
    return clamp(normalized, info["lower_bound"], info["upper_bound"])
end

"""
    prepare_input_tensor(jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                        jets::Vector{EEJet}, 
                        config::Dict, 
                        feature_data::Dict) -> Dict{String, Array}

Prepare input tensors for the neural network from jet constituents.

# Arguments
- `jcs`: Vector of jet constituents (structured as a vector of StructVector of ReconstructedParticle)
- `jets`: Vector of jets (EEJet)
- `config`: JSON configuration for preprocessing
- `feature_data`: Dictionary containing all extracted features

# Returns
Dictionary of input tensors
"""
function prepare_input_tensor(jcs::Vector{<:JetConstituents},
                              jets::Vector{EEJet},
                              config::Dict,
                              feature_data::Dict)

    # Get input names and variable info
    input_names = config["input_names"]

    # Initialize input tensor dictionary
    input_tensors = Dict{String, Array{Float32}}()

    # Get max length for padding
    max_length = config["pf_points"]["var_length"]

    # Process each jet
    for (i, jet) in enumerate(jets)
        # Initialize tensors
        if i == 1
            # Initialize tensor shapes based on config
            for input_name in input_names
                if input_name == "pf_features"
                    feature_vars = length(config[input_name]["var_names"])
                    input_tensors[input_name] = zeros(Float32, 1, feature_vars, max_length)
                    # elseif input_name == "pf_points"
                    #     points_vars = length(config[input_name]["var_names"])
                    #     input_tensors[input_name] = zeros(Float32, 1, points_vars, max_length)
                elseif input_name == "pf_vectors"
                    vector_vars = length(config[input_name]["var_names"])
                    input_tensors[input_name] = zeros(Float32, 1, vector_vars, max_length)
                elseif input_name == "pf_mask"
                    input_tensors[input_name] = zeros(Float32, 1, 1, max_length)
                end
            end
        end

        # Fill each tensor for this jet
        constituents = jcs[i]
        num_constituents = min(length(constituents), max_length)

        # Fill mask (1 for valid constituents, 0 for padding)
        if haskey(feature_data, "pf_mask")
            for j in 1:num_constituents
                input_tensors["pf_mask"][1, 1, j] = 1.0f0
            end
        end

        # Fill points
        if haskey(feature_data, "pf_points") && haskey(input_tensors, "pf_points")
            for (var_idx, var_name) in enumerate(config["pf_points"]["var_names"])
                var_info = config["pf_points"]["var_infos"][var_name]

                for j in 1:num_constituents
                    if j <= length(feature_data["pf_points"][var_name][i])
                        raw_value = feature_data["pf_points"][var_name][i][j]
                        norm_value = normalize_feature(raw_value, var_info)
                        input_tensors["pf_points"][1, var_idx, j] = norm_value
                    end
                end
            end
        end

        # Fill features
        if haskey(feature_data, "pf_features") && haskey(input_tensors, "pf_features")
            for (var_idx, var_name) in enumerate(config["pf_features"]["var_names"])
                var_info = config["pf_features"]["var_infos"][var_name]

                for j in 1:num_constituents
                    if haskey(feature_data["pf_features"], var_name) &&
                       j <= length(feature_data["pf_features"][var_name][i])
                        raw_value = feature_data["pf_features"][var_name][i][j]
                        norm_value = normalize_feature(raw_value, var_info)
                        input_tensors["pf_features"][1, var_idx, j] = norm_value
                    end
                end
            end
        end

        # Fill vectors (energies, momenta)
        if haskey(feature_data, "pf_vectors") && haskey(input_tensors, "pf_vectors")
            for (var_idx, var_name) in enumerate(config["pf_vectors"]["var_names"])
                for j in 1:num_constituents
                    if haskey(feature_data["pf_vectors"], var_name) &&
                       j <= length(feature_data["pf_vectors"][var_name][i])
                        input_tensors["pf_vectors"][1, var_idx, j] = feature_data["pf_vectors"][var_name][i][j]
                    end
                end
            end
        end
    end

    return input_tensors
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
Vector of flavour probabilities for each jet
"""
function get_weights(slot::Int, vars::Dict{String, Dict{String, Vector{Vector{Float32}}}},
                     jets::Vector{EEJet}, jcs::Vector{<:JetConstituents},
                     json_config::Dict, model::ONNXRunTime.InferenceSession)

    # Prepare input tensor
    input_tensors = prepare_input_tensor(jcs, jets, json_config, vars)

    # Run inference
    output = model(input_tensors)

    # Extract probabilities
    probabilities = output["softmax"]

    # Convert to desired output format (one vector per jet)
    result = Vector{Vector{Float32}}()

    # Reshape output to get probabilities for each jet
    num_classes = size(probabilities, 2)
    for i in 1:length(jets)
        jet_probs = Vector{Float32}(undef, num_classes)
        for c in 1:num_classes
            jet_probs[c] = probabilities[1, c]
        end
        push!(result, jet_probs)
    end

    return result
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
function get_weight(jet_weights::Vector{Vector{Float32}}, weight_idx::Int)
    if weight_idx < 0
        error("Invalid index requested for jet flavour weight.")
    end

    result = Vector{Float32}()

    for jet_weight in jet_weights
        if weight_idx >= length(jet_weight)
            error("Flavour weight index exceeds the number of weights registered.")
        end

        push!(result, jet_weight[weight_idx + 1])  # +1 for Julia's 1-based indexing
    end

    return result
end

"""
    inference(json_config_path::AbstractString, onnx_model_path::AbstractString, df::DataFrame,
                jets::Vector{EEJet}, jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                feature_data::Dict) -> DataFrame

Run flavour tagging inference on a collection of jets.

# Arguments
- `json_config_path`: Path to the JSON configuration file
- `onnx_model_path`: Path to the ONNX model file
- `jets`: Vector of jets
- `jcs`: Vector of jet constituents
- `feature_data`: Dictionary containing all extracted features

# Returns
DataFrame with added flavour tagging scores
"""
function inference(json_config_path::AbstractString, onnx_model_path::AbstractString,
                   jets::Vector{EEJet},
                   jcs::Vector{StructVector{EDM4hep.ReconstructedParticle}},
                   feature_data::Dict)

    # Extract input variables/score names from JSON file
    initvars = String[]
    variables = String[]
    scores = String[]

    config = JSON.parsefile(json_config_path)

    # Extract feature names
    for varname in config["pf_features"]["var_names"]
        push!(initvars, varname)
        push!(variables, varname)
    end

    # Extract vector names
    for varname in config["pf_vectors"]["var_names"]
        push!(initvars, varname)
        push!(variables, varname)
    end

    # Extract output names
    for scorename in config["output_names"]
        push!(scores, scorename)
    end

    # Setup model
    model, _ = setup_onnx_runtime(onnx_model_path, json_config_path, initvars)

    # Run inference
    weights = get_weights(0, feature_data, jets, jcs, config, model)

    # Extract individual scores
    jet_scores = Dict{String, Vector{Float32}}()

    for (i, scorename) in enumerate(scores)
        jet_scores[scorename] = get_weight(weights, i - 1)  # Adjust for 0-based indexing in get_weight
    end

    return jet_scores
end

# TODO: Add primary vertex as an argument (from MC Particle)
"""
    extract_features(jets::Vector{EEJet}, jcs::Vector{<:JetConstituents}, 
                    tracks::AbstractVector{EDM4hep.TrackState}, bz::Float32, 
                    track_L::AbstractArray{T} where T <: AbstractFloat, 
                    trackdata::AbstractVector{EDM4hep.Track}=AbstractVector{EDM4hep.Track}(), 
                    trackerhits::AbstractVector{EDM4hep.TrackerHit}=AbstractVector{EDM4hep.TrackerHit}(), 
                    gammadata::AbstractVector{EDM4hep.Cluster}=AbstractVector{EDM4hep.Cluster}(), 
                    nhdata::AbstractVector{EDM4hep.Cluster}=AbstractVector{EDM4hep.Cluster}(), 
                    calohits::AbstractVector{EDM4hep.CalorimeterHit}=AbstractVector{EDM4hep.CalorimeterHit}(), 
                    dNdx::AbstractVector{EDM4hep.Quantity}=AbstractVector{EDM4hep.Quantity}()) -> Dict

Extract all required features for jet flavour tagging.

# Arguments 
- `jets`: Vector of jets (EEJet)
- `jcs`: Vector of jet constituents
- `tracks`: StructVector of track states
- `bz`: Magnetic field strength
- `track_L`: Array of track lengths
- `trackdata`: Vector of track data (optional)
- `trackerhits`: Vector of tracker hits (optional)
- `gammadata`: Vector of gamma clusters (optional)
- `nhdata`: Vector of neutral hadron clusters (optional)
- `calohits`: Vector of calorimeter hits (optional)
- `dNdx`: Vector of dE/dx measurements (optional)   

# Returns
Dictionary containing all extracted features organized by input type
"""
function extract_features(jets::Vector{EEJet}, jcs::Vector{<:JetConstituents},
                          tracks::AbstractVector{EDM4hep.TrackState}, bz::Float32,
                          track_L::AbstractArray{T} where {T <: AbstractFloat},
                          trackdata::AbstractVector{EDM4hep.Track} = AbstractVector{EDM4hep.Track}(),
                          trackerhits::AbstractVector{EDM4hep.TrackerHit} = AbstractVector{EDM4hep.TrackerHit}(),
                          gammadata::AbstractVector{EDM4hep.Cluster} = AbstractVector{EDM4hep.Cluster}(),
                          nhdata::AbstractVector{EDM4hep.Cluster} = AbstractVector{EDM4hep.Cluster}(),
                          calohits::AbstractVector{EDM4hep.CalorimeterHit} = AbstractVector{EDM4hep.CalorimeterHit}(),
                          dNdx::AbstractVector{EDM4hep.Quantity} = AbstractVector{EDM4hep.Quantity}())

    # Primary vertex (0,0,0,0) for displacement calculations
    # TODO: Replace with actual primary vertex if available. Right now, the EDM4hep has bugs that don't allow me to get the f32 value.
    v_in = LorentzVector(0.0, 0.0, 0.0, 0.0)

    # Initialize feature containers
    features = Dict{String, Dict{String, Vector{Vector{Float32}}}}()

    # Initialize sub-dictionaries
    features["pf_points"] = Dict{String, Vector{Vector{Float32}}}()
    features["pf_features"] = Dict{String, Vector{Vector{Float32}}}()
    features["pf_vectors"] = Dict{String, Vector{Vector{Float32}}}()
    features["pf_mask"] = Dict{String, Vector{Vector{Float32}}}()

    # Extract basic features

    # Points (spatial coordinates)
    thetarel = get_thetarel_cluster(jets, jcs)
    phirel = get_phirel_cluster(jets, jcs)

    features["pf_points"]["pfcand_thetarel"] = thetarel
    features["pf_points"]["pfcand_phirel"] = phirel

    # Extract PF features
    erel_log = get_erel_log_cluster(jets, jcs)
    features["pf_features"]["pfcand_erel_log"] = erel_log
    features["pf_features"]["pfcand_thetarel"] = thetarel
    features["pf_features"]["pfcand_phirel"] = phirel

    # Track parameters and covariance matrices
    dxy = get_dxy(jcs, tracks, v_in, bz)
    dz = get_dz(jcs, tracks, v_in, bz)
    phi0 = get_phi0(jcs, tracks, v_in, bz)
    dptdpt = get_dptdpt(jcs, tracks)
    detadeta = get_detadeta(jcs, tracks)
    dphidphi = get_dphidphi(jcs, tracks)
    dxydxy = get_dxydxy(jcs, tracks)
    dzdz = get_dzdz(jcs, tracks)
    dxydz = get_dxydz(jcs, tracks)
    dphidxy = get_dphidxy(jcs, tracks)
    dlambdadz = get_dlambdadz(jcs, tracks)
    dxyc = get_dxyc(jcs, tracks)
    dxyctgtheta = get_dxyctgtheta(jcs, tracks)
    phic = get_phic(jcs, tracks)
    phidz = get_phidz(jcs, tracks)
    phictgtheta = get_phictgtheta(jcs, tracks)
    cdz = get_cdz(jcs, tracks)
    cctgtheta = get_cctgtheta(jcs, tracks)

    # Add track parameters to features
    features["pf_features"]["pfcand_dptdpt"] = dptdpt
    features["pf_features"]["pfcand_detadeta"] = detadeta
    features["pf_features"]["pfcand_dphidphi"] = dphidphi
    features["pf_features"]["pfcand_dxydxy"] = dxydxy
    features["pf_features"]["pfcand_dzdz"] = dzdz
    features["pf_features"]["pfcand_dxydz"] = dxydz
    features["pf_features"]["pfcand_dphidxy"] = dphidxy
    features["pf_features"]["pfcand_dlambdadz"] = dlambdadz
    features["pf_features"]["pfcand_dxyc"] = dxyc
    features["pf_features"]["pfcand_dxyctgtheta"] = dxyctgtheta
    features["pf_features"]["pfcand_phic"] = phic
    features["pf_features"]["pfcand_phidz"] = phidz
    features["pf_features"]["pfcand_phictgtheta"] = phictgtheta
    features["pf_features"]["pfcand_cdz"] = cdz
    features["pf_features"]["pfcand_cctgtheta"] = cctgtheta

    # Particle ID
    jcs_isChargedHad = get_is_charged_had(jcs)

    # Time-of-flight and dE/dx if data available
    if !isempty(track_L) && !isempty(trackdata) && !isempty(trackerhits) &&
       !isempty(gammadata) && !isempty(nhdata) && !isempty(calohits)
        mtof = get_mtof(jcs, track_L, trackdata, trackerhits, gammadata, nhdata, calohits,
                        v_in)
        features["pf_features"]["pfcand_mtof"] = mtof
    else
        # Empty vectors if data not available
        features["pf_features"]["pfcand_mtof"] = [Float32[] for _ in 1:length(jets)]
    end

    if !isempty(dNdx) && !isempty(trackdata)
        dndx_vals = get_dndx(jcs, dNdx, trackdata, jcs_isChargedHad)
        features["pf_features"]["pfcand_dndx"] = dndx_vals
    else
        features["pf_features"]["pfcand_dndx"] = [Float32[] for _ in 1:length(jets)]
    end

    # Particle type information
    charge = get_charge(jcs)
    isMu = get_is_mu(jcs)
    isEl = get_is_el(jcs)
    isChargedHad = jcs_isChargedHad
    isGamma = get_is_gamma(jcs)
    isNeutralHad = get_is_neutral_had(jcs)

    features["pf_features"]["pfcand_charge"] = charge
    features["pf_features"]["pfcand_isMu"] = isMu
    features["pf_features"]["pfcand_isEl"] = isEl
    features["pf_features"]["pfcand_isChargedHad"] = isChargedHad
    features["pf_features"]["pfcand_isGamma"] = isGamma
    features["pf_features"]["pfcand_isNeutralHad"] = isNeutralHad

    # Displacement variables
    features["pf_features"]["pfcand_dxy"] = dxy
    features["pf_features"]["pfcand_dz"] = dz

    # B-tagging variables
    btagSip2dVal = get_btagSip2dVal(jets, dxy, phi0, bz)
    btagSip2dSig = get_btagSip2dSig(btagSip2dVal, dxydxy)
    btagSip3dVal = get_btagSip3dVal(jets, dxy, dz, phi0, bz)
    btagSip3dSig = get_btagSip3dSig(btagSip3dVal, dxydxy, dzdz)
    btagJetDistVal = get_btagJetDistVal(jets, jcs, dxy, dz, phi0, bz)
    btagJetDistSig = get_btagJetDistSig(btagJetDistVal, dxydxy, dzdz)

    features["pf_features"]["pfcand_btagSip2dVal"] = btagSip2dVal
    features["pf_features"]["pfcand_btagSip2dSig"] = btagSip2dSig
    features["pf_features"]["pfcand_btagSip3dVal"] = btagSip3dVal
    features["pf_features"]["pfcand_btagSip3dSig"] = btagSip3dSig
    features["pf_features"]["pfcand_btagJetDistVal"] = btagJetDistVal
    features["pf_features"]["pfcand_btagJetDistSig"] = btagJetDistSig

    # Vector features (energy and momentum)
    e = get_e(jcs)
    p = get_p(jcs)

    features["pf_vectors"]["pfcand_e"] = e
    features["pf_vectors"]["pfcand_p"] = p

    # Add mask (all 1s for real particles, 0s for padding)
    mask = [fill(1.0f0, length(constituents)) for constituents in jcs]
    features["pf_mask"]["pfcand_mask"] = mask

    return features
end

end # module
