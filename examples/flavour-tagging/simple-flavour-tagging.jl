#!/usr/bin/env julia

"""
Simple Jet Flavour Tagging Example

This script demonstrates how to:
1. Load EDM4hep event data
2. Reconstruct jets using JetReconstruction
3. Extract features for flavour tagging
4. Run ONNX neural network inference
5. Get flavour probabilities for each jet

Run with: julia --project simple-flavour-tagging.jl
"""

using EDM4hep
using EDM4hep.RootIO
using LorentzVectorHEP
using JSON
using ONNXRunTime
using PhysicalConstants
using StructArrays
using JetReconstruction

function main()
    # Paths to model files
    model_dir = "data/wc_pt_7classes_12_04_2023"
    onnx_path = joinpath(model_dir, "fccee_flavtagging_edm4hep_wc_v1.onnx")
    json_path = joinpath(model_dir, "fccee_flavtagging_edm4hep_wc_v1.json")

    # Check if model files exist
    if !isfile(onnx_path)
        error("ONNX model not found at: $onnx_path")
    end
    if !isfile(json_path)
        error("JSON config not found at: $json_path")
    end

    println("Loading flavour tagging model...")
    config = JSON.parsefile(json_path)
    model = ONNXRunTime.load_inference(onnx_path)

    println("\nThe model predicts these flavour classes:")
    for class_name in config["output_names"]
        println("  - $class_name")
    end

    # Path to ROOT file with EDM4hep data
    edm4hep_path = "data/events_080263084.root"
    # edm4hep_path = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/wzp6_ee_nunuH_ecm240/events_080263084.root"
    if !isfile(edm4hep_path)
        error("EDM4hep data file not found at: $edm4hep_path")
    end

    println("\nLoading EDM4hep data...")
    reader = RootIO.Reader(edm4hep_path)
    events = RootIO.get(reader, "events")
    println("Loaded $(length(events)) events")

    # Process a specific event (event #12 as in the notebook)
    event_id = 15
    println("\nProcessing event #$event_id")
    evt = events[event_id]

    # Get reconstructed particles and tracks
    recps = RootIO.get(reader, evt, "ReconstructedParticles")
    tracks = RootIO.get(reader, evt, "EFlowTrack_1")

    # Get needed collections for feature extraction
    bz = RootIO.get(reader, evt, "magFieldBz", register = false)[1]
    trackdata = RootIO.get(reader, evt, "EFlowTrack")
    trackerhits = RootIO.get(reader, evt, "TrackerHits")
    gammadata = RootIO.get(reader, evt, "EFlowPhoton")
    nhdata = RootIO.get(reader, evt, "EFlowNeutralHadron")
    calohits = RootIO.get(reader, evt, "CalorimeterHits")
    dNdx = RootIO.get(reader, evt, "EFlowTrack_2")
    track_L = RootIO.get(reader, evt, "EFlowTrack_L", register = false)

    println("  - $(length(recps)) reconstructed particles")
    println("  - $(length(tracks)) tracks")
    println("  - Magnetic field Bz = $bz T")

    # Reconstruct jets
    println("\nReconstructing jets...")
    cs = jet_reconstruct(recps; p = 1.0, R = 2.0, algorithm = JetAlgorithm.EEKt)

    # Get 2 exclusive jets
    jets = exclusive_jets(cs; njets = 2, T = EEJet)
    println("Found $(length(jets)) jets")

    # Print jet properties
    for (i, jet) in enumerate(jets)
        println("\nJet $i:")
        println("  - Energy: $(round(jet.E, digits=2)) GeV")
        println("  - Pt: $(round(JetReconstruction.pt(jet), digits=2)) GeV")
        println("  - Eta: $(round(JetReconstruction.eta(jet), digits=3))")
        println("  - Phi: $(round(JetReconstruction.phi(jet), digits=3))")
        println("  - Mass: $(round(JetReconstruction.mass(jet), digits=2)) GeV")
    end

    # Get jet constituents
    println("\nExtracting jet constituents...")
    constituent_indices = [constituent_indexes(jet, cs) for jet in jets]

    # Access the extension module
    ext_mod = Base.get_extension(JetReconstruction, :JetFlavourTagging)
    if isnothing(ext_mod)
        error("JetFlavourTagging extension not loaded")
    end

    jet_constituents = ext_mod.build_constituents_cluster(recps, constituent_indices)

    for (i, constituents) in enumerate(jet_constituents)
        println("  - Jet $i has $(length(constituents)) constituents")
    end

    # Extract features for flavour tagging
    println("\nExtracting features for flavour tagging...")
    feature_data = JetReconstruction.extract_features(jets,
                                                      jet_constituents,
                                                      tracks,
                                                      bz,
                                                      track_L,
                                                      config,
                                                      trackdata,
                                                      trackerhits,
                                                      gammadata,
                                                      nhdata,
                                                      calohits,
                                                      dNdx)

    # Prepare input tensors
    println("Preparing input tensors...")
    input_tensors = JetReconstruction.prepare_input_tensor(jet_constituents,
                                                           jets,
                                                           config,
                                                           feature_data)

    # Run inference
    println("Running neural network inference...")
    weights = JetReconstruction.get_weights(0,  # Thread slot
                                            feature_data,
                                            jets,
                                            jet_constituents,
                                            config,
                                            model)

    # Extract and display results
    println("\n" * "="^60)
    println("FLAVOUR TAGGING RESULTS")
    println("="^60)

    for (jet_idx, jet) in enumerate(jets)
        println("\nJet $jet_idx (E=$(round(jet.E, digits=1)) GeV, pT=$(round(JetReconstruction.pt(jet), digits=1)) GeV):")
        println("-"^40)

        # Collect scores for this jet
        scores = Float32[]
        labels = String[]

        for (i, score_name) in enumerate(config["output_names"])
            score = JetReconstruction.get_weight(weights, i - 1)[jet_idx]
            push!(scores, score)
            push!(labels, score_name)
        end

        # Sort by probability (descending)
        sorted_indices = sortperm(scores, rev = true)

        # Display scores
        for idx in sorted_indices
            label = labels[idx]
            score = scores[idx]

            # Handle NaN or invalid scores
            if isnan(score) || isinf(score)
                flavor_map = Dict("recojet_isG" => "Gluon   ",
                                  "recojet_isQ" => "Light q ",
                                  "recojet_isS" => "Strange ",
                                  "recojet_isC" => "Charm   ",
                                  "recojet_isB" => "Bottom  ")
                formatted_label = get(flavor_map, label, label)
                println("  $formatted_label: [Invalid score]")
                continue
            end

            bar_length = Int(round(score * 30))
            bar = "█"^bar_length
            percentage = round(score * 100, digits = 1)

            # Format label
            flavor_map = Dict("recojet_isG" => "Gluon   ",
                              "recojet_isQ" => "Light q ",
                              "recojet_isS" => "Strange ",
                              "recojet_isC" => "Charm   ",
                              "recojet_isB" => "Bottom  ")

            formatted_label = get(flavor_map, label, label)
            println("  $formatted_label: $bar $(percentage)%")
        end

        # Identify most likely flavour
        max_idx = argmax(scores)
        max_label = labels[max_idx]
        max_score = scores[max_idx]

        flavour_name = Dict("recojet_isG" => "gluon",
                            "recojet_isQ" => "light quark",
                            "recojet_isS" => "strange",
                            "recojet_isC" => "charm",
                            "recojet_isB" => "bottom")[max_label]

        println("\n  → Most likely: $(flavour_name) ($(round(max_score * 100, digits=1))% confidence)")
    end

    println("\n" * "="^60)
    println("Processing complete!")
end

# Run the main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
