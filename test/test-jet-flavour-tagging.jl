# Tests for Jet Flavour Tagging Extension

include("common.jl")

# Check if the extension is available
const HAS_FLAVOUR_TAGGING = try
    using EDM4hep
    using ONNXRunTime
    using JSON
    true
catch
    false
end

if !HAS_FLAVOUR_TAGGING
    @warn "Jet Flavour Tagging extension dependencies not available, skipping tests"
else
    using EDM4hep
    using EDM4hep.RootIO
    using StaticArrays
    using LorentzVectorHEP
    using JetReconstruction
    using LoopVectorization
    using JSON
    using ONNXRunTime
    using StructArrays: StructVector

    @testset "Jet Flavour Tagging Extension" begin
        
        @testset "Extension Functions Defined" begin
            # Verify all exported functions are defined
            @test isdefined(JetReconstruction, :build_constituents_cluster)
            @test isdefined(JetReconstruction, :extract_features)
            @test isdefined(JetReconstruction, :setup_onnx_runtime)
            @test isdefined(JetReconstruction, :prepare_input_tensor)
            @test isdefined(JetReconstruction, :get_weights)
            @test isdefined(JetReconstruction, :get_weight)
        end

        # Navigate to the examples folder from the test folder
        package_root = dirname(@__DIR__)  # Go up from test/ to JetReconstruction.jl/
        onnx_data_dir = joinpath(package_root, "examples", "flavour-tagging", "data", "wc_pt_7classes_12_04_2023")
        edm4hep_dir = joinpath(package_root, "examples", "flavour-tagging", "data")
        
        # Paths to the actual files
        onnx_path = joinpath(onnx_data_dir, "fccee_flavtagging_edm4hep_wc_v1.onnx")
        json_path = joinpath(onnx_data_dir, "fccee_flavtagging_edm4hep_wc_v1.json")
        edm4hep_example_file = joinpath(edm4hep_dir, "events_080263084.root")

        # Check if files exist
        has_onnx_file = isfile(onnx_path)
        has_json_file = isfile(json_path)
        has_edm4hep_file = isfile(edm4hep_example_file)
        has_test_files = has_onnx_file && has_json_file && has_edm4hep_file
        
        if !has_test_files
            @warn "Test files not found in examples folder. Missing files:" *
                  (!has_onnx_file ? "\n  - $onnx_path" : "") *
                  (!has_json_file ? "\n  - $json_path" : "") *
                  (!has_edm4hep_file ? "\n  - $edm4hep_example_file" : "")
        end

        if has_test_files
            @testset "Full Flavour Tagging Pipeline" begin
                # Load the EDM4hep data
                reader = RootIO.Reader(edm4hep_example_file)
                events = RootIO.get(reader, "events")
                @test length(events) > 0
                
                # Use event #16 as in the notebook
                event_id = 16
                evt = events[event_id]
                
                # Get reconstructed particles and tracks
                recps = RootIO.get(reader, evt, "ReconstructedParticles")
                @test length(recps) > 0
                
                tracks = RootIO.get(reader, evt, "EFlowTrack_1")
                @test length(tracks) >= 0
                
                # Get needed collections for feature extraction
                bz = RootIO.get(reader, evt, "magFieldBz", register = false)[1]
                @test isa(bz, Float32)
                
                trackdata = RootIO.get(reader, evt, "EFlowTrack")
                trackerhits = RootIO.get(reader, evt, "TrackerHits")
                gammadata = RootIO.get(reader, evt, "EFlowPhoton")
                nhdata = RootIO.get(reader, evt, "EFlowNeutralHadron")
                calohits = RootIO.get(reader, evt, "CalorimeterHits")
                dNdx = RootIO.get(reader, evt, "EFlowTrack_2")
                track_L = RootIO.get(reader, evt, "EFlowTrack_L", register = false)
                
                @testset "Jet Reconstruction" begin
                    # Cluster jets using the EEkt algorithm
                    cs = jet_reconstruct(recps; p = 1.0, R = 2.0, algorithm = JetAlgorithm.EEKt)
                    @test isa(cs, ClusterSequence)
                    
                    # Get 2 exclusive jets
                    jets = exclusive_jets(cs; njets=2, T=EEJet)
                    @test length(jets) == 2
                    @test isa(jets[1], EEJet)
                end
                
                @testset "build_constituents_cluster" begin
                    cs = jet_reconstruct(recps; p = 1.0, R = 2.0, algorithm = JetAlgorithm.EEKt)
                    jets = exclusive_jets(cs; njets=2, T=EEJet)
                    
                    # Get constituent indices for each jet
                    constituent_indices = [constituent_indexes(jet, cs) for jet in jets]
                    @test length(constituent_indices) == 2
                    @test all(indices -> length(indices) > 0, constituent_indices)
                    
                    # Build constituents
                    jet_constituents = JetReconstruction.build_constituents_cluster(recps, constituent_indices)
                    @test length(jet_constituents) == 2
                    @test all(jc -> isa(jc, StructVector{ReconstructedParticle}), jet_constituents)
                end
                
                @testset "extract_features" begin
                    cs = jet_reconstruct(recps; p = 1.0, R = 2.0, algorithm = JetAlgorithm.EEKt)
                    jets = exclusive_jets(cs; njets=2, T=EEJet)
                    constituent_indices = [constituent_indexes(jet, cs) for jet in jets]
                    jet_constituents = JetReconstruction.build_constituents_cluster(recps, constituent_indices)
                    
                    # Extract features
                    feature_data = JetReconstruction.extract_features(
                        jets, 
                        jet_constituents, 
                        tracks, 
                        bz, 
                        track_L, 
                        trackdata, 
                        trackerhits, 
                        gammadata, 
                        nhdata, 
                        calohits, 
                        dNdx
                    )
                    
                    @test isa(feature_data, Dict)
                    @test haskey(feature_data, "Particles")
                    @test haskey(feature_data, "Tracks")
                end
                
                @testset "setup_onnx_runtime" begin
                    # Test ONNX runtime setup
                    model = JetReconstruction.setup_onnx_runtime(onnx_path, json_path)
                    @test isa(model, ONNXRunTime.InferenceSession)
                    
                    # Also test loading config separately
                    config = JSON.parsefile(json_path)
                    @test isa(config, Dict)
                    @test haskey(config, "output_names")
                    @test length(config["output_names"]) == 5  # 5 flavor classes
                end
                
                @testset "prepare_input_tensor" begin
                    # Setup everything needed
                    cs = jet_reconstruct(recps; p = 1.0, R = 2.0, algorithm = JetAlgorithm.EEKt)
                    jets = exclusive_jets(cs; njets=2, T=EEJet)
                    constituent_indices = [constituent_indexes(jet, cs) for jet in jets]
                    jet_constituents = JetReconstruction.build_constituents_cluster(recps, constituent_indices)
                    
                    feature_data = JetReconstruction.extract_features(
                        jets, jet_constituents, tracks, bz, track_L,
                        trackdata, trackerhits, gammadata, nhdata, calohits, dNdx
                    )
                    
                    config = JSON.parsefile(json_path)
                    
                    # Prepare input tensors
                    input_tensors = JetReconstruction.prepare_input_tensor(
                        jet_constituents, jets, config, feature_data
                    )
                    
                    @test isa(input_tensors, Dict)
                    @test all(v -> isa(v, Array), values(input_tensors))
                end
                
                @testset "get_weights and get_weight" begin
                    # Full pipeline to get weights
                    cs = jet_reconstruct(recps; p = 1.0, R = 2.0, algorithm = JetAlgorithm.EEKt)
                    jets = exclusive_jets(cs; njets=2, T=EEJet)
                    constituent_indices = [constituent_indexes(jet, cs) for jet in jets]
                    jet_constituents = JetReconstruction.build_constituents_cluster(recps, constituent_indices)
                    
                    feature_data = JetReconstruction.extract_features(
                        jets, jet_constituents, tracks, bz, track_L,
                        trackdata, trackerhits, gammadata, nhdata, calohits, dNdx
                    )
                    
                    config = JSON.parsefile(json_path)
                    model = JetReconstruction.setup_onnx_runtime(onnx_path, json_path)
                    
                    # Get weights
                    weights = JetReconstruction.get_weights(
                        0,  # Thread slot
                        feature_data,
                        jets,
                        jet_constituents,
                        config,
                        model
                    )
                    
                    @test isa(weights, Vector{Vector{Float32}})
                    @test length(weights) == 2  # One weight vector per jet
                    @test all(w -> length(w) == 5, weights)  # 5 flavor classes per jet
                    
                    # Test get_weight for extracting individual scores
                    for i in 0:4  # 5 flavor classes (0-indexed)
                        score = JetReconstruction.get_weight(weights, i)
                        @test isa(score, Vector{Float32})
                        @test length(score) == 2  # One score per jet
                        @test all(s -> 0.0 <= s <= 1.0, score)  # Probabilities should be in [0,1]
                    end
                    
                    # Test that probabilities sum to approximately 1
                    for jet_idx in 1:2
                        prob_sum = sum(weights[jet_idx])
                        @test isapprox(prob_sum, 1.0, atol=1e-4)
                    end
                end
            end
        end
    end
end