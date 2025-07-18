module JetFlavourHelper

using JSON
using ONNXRunTime
using StructArrays: StructVector
using JetReconstruction
using EDM4hep
using LorentzVectorHEP

# Utility functions for the jet features from EDM4hep
include("JetConstituentUtils.jl")

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
    prepare_input_tensor(jets_constituents::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                        jets::Vector{EEJet}, 
                        config::Dict, 
                        feature_data::Dict,
                        jet_index::Int=1) -> Dict{String, Array}

Prepare input tensors for the neural network from jet constituents.

# Arguments
- `jets_constituents`: Vector of jet constituents (structured as a vector of StructVector of ReconstructedParticle)
- `jets`: Vector of jets (EEJet)
- `config`: JSON configuration for preprocessing
- `feature_data`: Dictionary containing all extracted features
- `jet_index`: Index of the jet to process (default 1)

# Returns
Dictionary of input tensors
"""
function prepare_input_tensor(jets_constituents::Vector{<:JetConstituents},
                              jets::Vector{EEJet},
                              config::Dict,
                              feature_data::Dict,
                              jet_index::Int = 1)

    # Get input names and variable info
    input_names = config["input_names"]

    # Initialize input tensor dictionary
    input_tensors = Dict{String, Array{Float32}}()

    # Get max length for padding
    max_length = config["pf_points"]["var_length"]

    # Initialize tensors for single jet processing
    for input_name in input_names
        if input_name == "pf_features"
            feature_vars = length(config[input_name]["var_names"])
            input_tensors[input_name] = zeros(Float32, 1, feature_vars, max_length)
        elseif input_name == "pf_vectors"
            vector_vars = length(config[input_name]["var_names"])
            input_tensors[input_name] = zeros(Float32, 1, vector_vars, max_length)
        elseif input_name == "pf_mask"
            input_tensors[input_name] = zeros(Float32, 1, 1, max_length)
        end
    end

    # Note: This function prepares tensors for a single jet at a time
    # The caller can loop through jets and process them individually
    # A batch processing version is planned for next version under the inference() function.

    # Process the specified jet
    if length(jets) >= jet_index && jet_index > 0
        i = jet_index
        jet = jets[1]  # We're processing single jets, so always use first element

        # Fill each tensor for this jet
        constituents = jets_constituents[1]  # Single jet passed
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
                     jets::Vector{EEJet}, jets_constituents::Vector{<:JetConstituents},
                     json_config::Dict, model::ONNXRunTime.InferenceSession)

    # The model processes one jet at a time
    result = Vector{Vector{Float32}}()

    # Process each jet individually
    for i in 1:length(jets)
        # Create single-jet arrays
        single_jet = [jets[i]]
        single_constituents = [jets_constituents[i]]

        # Create single-jet feature data by extracting only features for this jet
        single_jet_vars = Dict{String, Dict{String, Vector{Vector{Float32}}}}()
        for (category, features) in vars
            single_jet_vars[category] = Dict{String, Vector{Vector{Float32}}}()
            for (fname, fvalues) in features
                # Extract only the features for jet i
                if i <= length(fvalues)
                    single_jet_vars[category][fname] = [fvalues[i]]
                else
                    # If no features for this jet, create empty array
                    single_jet_vars[category][fname] = [Float32[]]
                end
            end
        end

        # Prepare input tensor for this single jet with extracted features
        input_tensors = prepare_input_tensor(single_constituents, single_jet, json_config,
                                             single_jet_vars, 1)

        # Run inference
        output = model(input_tensors)

        # Extract probabilities
        probabilities = output["softmax"]

        # Get probabilities for this jet
        num_classes = size(probabilities, 2)
        jet_probs = Vector{Float32}(undef, num_classes)
        for c in 1:num_classes
            jet_probs[c] = probabilities[1, c]  # Always index 1 since we process one jet at a time
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
                jets::Vector{EEJet}, jets_constituents::Vector{StructVector{EDM4hep.ReconstructedParticle}}, 
                feature_data::Dict) -> DataFrame

Run flavour tagging inference on a collection of jets.

# Arguments
- `json_config_path`: Path to the JSON configuration file
- `onnx_model_path`: Path to the ONNX model file
- `jets`: Vector of jets
- `jets_constituents`: Vector of jet constituents
- `feature_data`: Dictionary containing all extracted features

# Returns
DataFrame with added flavour tagging scores
"""
function inference(json_config_path::AbstractString, onnx_model_path::AbstractString,
                   jets::Vector{EEJet},
                   jets_constituents::Vector{StructVector{EDM4hep.ReconstructedParticle}},
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
    weights = get_weights(0, feature_data, jets, jets_constituents, config, model)

    # Extract individual scores
    jet_scores = Dict{String, Vector{Float32}}()

    for (i, scorename) in enumerate(scores)
        jet_scores[scorename] = get_weight(weights, i - 1)  # Adjust for 0-based indexing in get_weight
    end

    return jet_scores
end

# TODO: Add primary vertex as an argument (from MC Particle)
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
function extract_features(jets::Vector{EEJet}, jets_constituents::Vector{<:JetConstituents},
                          tracks::AbstractVector{EDM4hep.TrackState}, bz::Float32,
                          track_L::AbstractArray{T} where {T <: AbstractFloat},
                          json_config::Dict,
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

    # Cache for computed features to avoid redundant calculations
    computed_features = Dict{String, Vector{Vector{Float32}}}()

    # Create kwargs dict with all available data
    kwargs = Dict(:jets => jets,
                  :tracks => tracks,
                  :v_in => v_in,
                  :bz => bz,
                  :track_L => track_L,
                  :trackdata => trackdata,
                  :trackerhits => trackerhits,
                  :gammadata => gammadata,
                  :nhdata => nhdata,
                  :calohits => calohits,
                  :dNdx => dNdx)

    # Process each input type specified in the JSON
    for input_name in json_config["input_names"]
        if !haskey(json_config, input_name) || !haskey(json_config[input_name], "var_names")
            continue
        end

        features[input_name] = Dict{String, Vector{Vector{Float32}}}()

        # Extract each requested variable
        for var_name in json_config[input_name]["var_names"]
            # Remove "pfcand_" prefix if present
            clean_name = replace(var_name, "pfcand_" => "")

            # Check if already computed
            if haskey(computed_features, clean_name)
                features[input_name][var_name] = computed_features[clean_name]
                continue
            end

            # Get the extractor function
            extractor = get_feature_extractor(clean_name)
            if extractor === nothing
                # Feature not implemented, create empty arrays
                features[input_name][var_name] = [Float32[] for _ in 1:length(jets)]
                continue
            end

            # Handle dependencies for certain features
            deps = get_feature_dependencies(clean_name)
            if !isempty(deps)
                dep_kwargs = Dict{Symbol, Any}()

                for (dep_key, dep_name) in deps
                    if !haskey(computed_features, dep_name)
                        # Recursively compute dependency
                        dep_extractor = get_feature_extractor(dep_name)
                        if dep_extractor !== nothing
                            # Check if this dependency has its own dependencies
                            sub_deps = get_feature_dependencies(dep_name)
                            if !isempty(sub_deps)
                                sub_dep_kwargs = Dict{Symbol, Any}()
                                for (sub_dep_key, sub_dep_name) in sub_deps
                                    if !haskey(computed_features, sub_dep_name)
                                        sub_dep_extractor = get_feature_extractor(sub_dep_name)
                                        if sub_dep_extractor !== nothing
                                            computed_features[sub_dep_name] = sub_dep_extractor(jets_constituents;
                                                                                                kwargs...)
                                        end
                                    end
                                    if haskey(computed_features, sub_dep_name)
                                        sub_dep_kwargs[sub_dep_key] = computed_features[sub_dep_name]
                                    end
                                end
                                # Compute dependency with its dependencies
                                temp_kwargs = merge(kwargs, sub_dep_kwargs)
                                computed_features[dep_name] = dep_extractor(jets_constituents;
                                                                            temp_kwargs...)
                            else
                                # No sub-dependencies, compute directly
                                computed_features[dep_name] = dep_extractor(jets_constituents;
                                                                            kwargs...)
                            end
                        end
                    end
                    if haskey(computed_features, dep_name)
                        dep_kwargs[dep_key] = computed_features[dep_name]
                    end
                end

                # Merge dependency kwargs with main kwargs
                merge!(kwargs, dep_kwargs)
            end

            # Special handling for features that need isChargedHad
            if clean_name == "dndx" && !haskey(kwargs, :jets_constituents_isChargedHad)
                if !haskey(computed_features, "isChargedHad")
                    computed_features["isChargedHad"] = get_is_charged_had(jets_constituents)
                end
                kwargs[:jets_constituents_isChargedHad] = computed_features["isChargedHad"]
            end

            # Extract the feature
            try
                feature_data = extractor(jets_constituents; kwargs...)
                computed_features[clean_name] = feature_data
                features[input_name][var_name] = feature_data
            catch e
                # If extraction fails, create empty arrays
                @warn "Failed to extract feature $var_name: $e"
                features[input_name][var_name] = [Float32[] for _ in 1:length(jets)]
            end
        end
    end

    return features
end

"""
    get_feature_extractor(feature_name::String)

Return the appropriate feature extraction function for a given feature name.

# Arguments
- `feature_name`: Name of the feature (without "pfcand_" prefix)

# Returns
A function that can extract the specified feature, or nothing if not found
"""
function get_feature_extractor(feature_name::String)
    # Map feature names to extraction functions
    feature_map = Dict(
                       # Basic kinematic features
                       "e" => (jets_constituents; kwargs...) -> get_e(jets_constituents),
                       "p" => (jets_constituents; kwargs...) -> get_p(jets_constituents),
                       "theta" => (jets_constituents; kwargs...) -> get_theta(jets_constituents),
                       "phi" => (jets_constituents; kwargs...) -> get_phi(jets_constituents),
                       "charge" => (jets_constituents; kwargs...) -> get_charge(jets_constituents),
                       "type" => (jets_constituents; kwargs...) -> get_type(jets_constituents),

                       # Relative kinematic features
                       "erel_log" => (jets_constituents; jets = nothing, kwargs...) -> get_erel_log_cluster(jets,
                                                                                                            jets_constituents),
                       "thetarel" => (jets_constituents; jets = nothing, kwargs...) -> get_thetarel_cluster(jets,
                                                                                                            jets_constituents),
                       "phirel" => (jets_constituents; jets = nothing, kwargs...) -> get_phirel_cluster(jets,
                                                                                                        jets_constituents),

                       # Particle ID features
                       "isMu" => (jets_constituents; kwargs...) -> get_is_mu(jets_constituents),
                       "isEl" => (jets_constituents; kwargs...) -> get_is_el(jets_constituents),
                       "isChargedHad" => (jets_constituents; kwargs...) -> get_is_charged_had(jets_constituents),
                       "isGamma" => (jets_constituents; kwargs...) -> get_is_gamma(jets_constituents),
                       "isNeutralHad" => (jets_constituents; kwargs...) -> get_is_neutral_had(jets_constituents),

                       # Track parameters
                       "dxy" => (jets_constituents; tracks = nothing, v_in = nothing, bz = nothing, kwargs...) -> get_dxy(jets_constituents,
                                                                                                                          tracks,
                                                                                                                          v_in,
                                                                                                                          bz),
                       "dz" => (jets_constituents; tracks = nothing, v_in = nothing, bz = nothing, kwargs...) -> get_dz(jets_constituents,
                                                                                                                        tracks,
                                                                                                                        v_in,
                                                                                                                        bz),
                       "phi0" => (jets_constituents; tracks = nothing, v_in = nothing, bz = nothing, kwargs...) -> get_phi0(jets_constituents,
                                                                                                                            tracks,
                                                                                                                            v_in,
                                                                                                                            bz),

                       # Covariance matrix elements
                       "dptdpt" => (jets_constituents; tracks = nothing, kwargs...) -> get_dptdpt(jets_constituents,
                                                                                                  tracks),
                       "detadeta" => (jets_constituents; tracks = nothing, kwargs...) -> get_detadeta(jets_constituents,
                                                                                                      tracks),
                       "dphidphi" => (jets_constituents; tracks = nothing, kwargs...) -> get_dphidphi(jets_constituents,
                                                                                                      tracks),
                       "dxydxy" => (jets_constituents; tracks = nothing, kwargs...) -> get_dxydxy(jets_constituents,
                                                                                                  tracks),
                       "dzdz" => (jets_constituents; tracks = nothing, kwargs...) -> get_dzdz(jets_constituents,
                                                                                              tracks),
                       "dxydz" => (jets_constituents; tracks = nothing, kwargs...) -> get_dxydz(jets_constituents,
                                                                                                tracks),
                       "dphidxy" => (jets_constituents; tracks = nothing, kwargs...) -> get_dphidxy(jets_constituents,
                                                                                                    tracks),
                       "dlambdadz" => (jets_constituents; tracks = nothing, kwargs...) -> get_dlambdadz(jets_constituents,
                                                                                                        tracks),
                       "dxyc" => (jets_constituents; tracks = nothing, kwargs...) -> get_dxyc(jets_constituents,
                                                                                              tracks),
                       "dxyctgtheta" => (jets_constituents; tracks = nothing, kwargs...) -> get_dxyctgtheta(jets_constituents,
                                                                                                            tracks),
                       "phic" => (jets_constituents; tracks = nothing, kwargs...) -> get_phic(jets_constituents,
                                                                                              tracks),
                       "phidz" => (jets_constituents; tracks = nothing, kwargs...) -> get_phidz(jets_constituents,
                                                                                                tracks),
                       "phictgtheta" => (jets_constituents; tracks = nothing, kwargs...) -> get_phictgtheta(jets_constituents,
                                                                                                            tracks),
                       "cdz" => (jets_constituents; tracks = nothing, kwargs...) -> get_cdz(jets_constituents,
                                                                                            tracks),
                       "cctgtheta" => (jets_constituents; tracks = nothing, kwargs...) -> get_cctgtheta(jets_constituents,
                                                                                                        tracks),

                       # B-tagging features
                       "btagSip2dVal" => (jets_constituents; jets = nothing, dxy = nothing, phi0 = nothing, bz = nothing, kwargs...) -> get_btagSip2dVal(jets,
                                                                                                                                                         dxy,
                                                                                                                                                         phi0,
                                                                                                                                                         bz),
                       "btagSip2dSig" => (jets_constituents; btagSip2dVal = nothing, dxydxy = nothing, kwargs...) -> get_btagSip2dSig(btagSip2dVal,
                                                                                                                                      dxydxy),
                       "btagSip3dVal" => (jets_constituents; jets = nothing, dxy = nothing, dz = nothing, phi0 = nothing, bz = nothing, kwargs...) -> get_btagSip3dVal(jets,
                                                                                                                                                                       dxy,
                                                                                                                                                                       dz,
                                                                                                                                                                       phi0,
                                                                                                                                                                       bz),
                       "btagSip3dSig" => (jets_constituents; btagSip3dVal = nothing, dxydxy = nothing, dzdz = nothing, kwargs...) -> get_btagSip3dSig(btagSip3dVal,
                                                                                                                                                      dxydxy,
                                                                                                                                                      dzdz),
                       "btagJetDistVal" => (jets_constituents; jets = nothing, dxy = nothing, dz = nothing, phi0 = nothing, bz = nothing, kwargs...) -> get_btagJetDistVal(jets,
                                                                                                                                                                           jets_constituents,
                                                                                                                                                                           dxy,
                                                                                                                                                                           dz,
                                                                                                                                                                           phi0,
                                                                                                                                                                           bz),
                       "btagJetDistSig" => (jets_constituents; btagJetDistVal = nothing, dxydxy = nothing, dzdz = nothing, kwargs...) -> get_btagJetDistSig(btagJetDistVal,
                                                                                                                                                            dxydxy,
                                                                                                                                                            dzdz),

                       # Special features
                       "mtof" => (jets_constituents; track_L = nothing, trackdata = nothing, trackerhits = nothing, gammadata = nothing, nhdata = nothing, calohits = nothing, v_in = nothing, kwargs...) -> get_mtof(jets_constituents,
                                                                                                                                                                                                                      track_L,
                                                                                                                                                                                                                      trackdata,
                                                                                                                                                                                                                      trackerhits,
                                                                                                                                                                                                                      gammadata,
                                                                                                                                                                                                                      nhdata,
                                                                                                                                                                                                                      calohits,
                                                                                                                                                                                                                      v_in),
                       "dndx" => (jets_constituents; dNdx = nothing, trackdata = nothing, jets_constituents_isChargedHad = nothing, kwargs...) -> get_dndx(jets_constituents,
                                                                                                                                                           dNdx,
                                                                                                                                                           trackdata,
                                                                                                                                                           jets_constituents_isChargedHad),

                       # Mask
                       "mask" => (jets_constituents; kwargs...) -> [fill(1.0f0,
                                                                         length(constituents))
                                                                    for constituents in jets_constituents])

    return get(feature_map, feature_name, nothing)
end

"""
    get_feature_dependencies(feature_name::String) -> Dict{Symbol, String}

Get the dependencies for features that require other computed features.

# Arguments
- `feature_name`: Name of the feature

# Returns
Dictionary mapping parameter names to feature names
"""
function get_feature_dependencies(feature_name::String)
    deps = Dict{Symbol, String}()

    if feature_name == "btagSip2dSig"
        deps[:btagSip2dVal] = "btagSip2dVal"
        deps[:dxydxy] = "dxydxy"
    elseif feature_name == "btagSip3dSig"
        deps[:btagSip3dVal] = "btagSip3dVal"
        deps[:dxydxy] = "dxydxy"
        deps[:dzdz] = "dzdz"
    elseif feature_name == "btagJetDistSig"
        deps[:btagJetDistVal] = "btagJetDistVal"
        deps[:dxydxy] = "dxydxy"
        deps[:dzdz] = "dzdz"
    elseif feature_name == "btagSip2dVal"
        deps[:dxy] = "dxy"
        deps[:phi0] = "phi0"
    elseif feature_name == "btagSip3dVal"
        deps[:dxy] = "dxy"
        deps[:dz] = "dz"
        deps[:phi0] = "phi0"
    elseif feature_name == "btagJetDistVal"
        deps[:dxy] = "dxy"
        deps[:dz] = "dz"
        deps[:phi0] = "phi0"
    end

    return deps
end

end # module
