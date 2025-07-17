# Tests for Jet Flavour Tagging Constituent Utilities

include("common.jl")

# Check if the extension is available
const HAS_FLAVOUR_TAGGING = try
    using EDM4hep
    using ONNXRunTime
    using PhysicalConstants
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
    using LorentzVectorHEP
    using ONNXRunTime
    using PhysicalConstants
    using JSON
    using StructArrays
    using JetReconstruction

    @testset "Jet Flavour Constituent Utilities" begin
        # Navigate to find the test reference file
        reference_file = joinpath(@__DIR__, "data", "jet-flavour-utils-event15.json")
        
        # Check if reference file exists
        if !isfile(reference_file)
            @warn "Reference output file not found at: $reference_file" *
                  "\nRun generate-jet-constituent-utils-outputs.jl first to generate it."
            @test_skip "Reference outputs not available"
        else
            # Load reference outputs
            reference_outputs = JSON.parsefile(reference_file)
            
            # Load the same data used to generate reference
            package_root = dirname(@__DIR__)
            edm4hep_dir = joinpath(package_root, "examples", "flavour-tagging", "data")
            edm4hep_example_file = joinpath(edm4hep_dir, "events_080263084.root")
            
            if !isfile(edm4hep_example_file)
                @test_skip "EDM4hep example file not found"
            else
                # Load data
                reader = RootIO.Reader(edm4hep_example_file)
                events = RootIO.get(reader, "events")
                event_id = reference_outputs["input_info"]["event_id"]
                evt = events[event_id]
                
                # Get all necessary data
                recps = RootIO.get(reader, evt, "ReconstructedParticles")
                tracks = RootIO.get(reader, evt, "EFlowTrack_1")
                bz = Float32(reference_outputs["input_info"]["bz"])
                trackdata = RootIO.get(reader, evt, "EFlowTrack")
                trackerhits = RootIO.get(reader, evt, "TrackerHits")
                gammadata = RootIO.get(reader, evt, "EFlowPhoton")
                nhdata = RootIO.get(reader, evt, "EFlowNeutralHadron")
                calohits = RootIO.get(reader, evt, "CalorimeterHits")
                dNdx = RootIO.get(reader, evt, "EFlowTrack_2")
                track_L = RootIO.get(reader, evt, "EFlowTrack_L", register = false)
                
                # Create vertex from reference
                v_info = reference_outputs["input_info"]["vertex"]
                V = LorentzVector(v_info["x"], v_info["y"], v_info["z"], v_info["t"])
                
                # Access the extension module first
                ext_mod = Base.get_extension(JetReconstruction, :JetFlavourTagging)
                @test !isnothing(ext_mod)
                
                # Access the helper module within the extension
                helper_mod = ext_mod.JetFlavourHelper
                
                # Reconstruct jets
                cs = jet_reconstruct(recps; p = 1.0, R = 2.0, algorithm = JetAlgorithm.EEKt)
                jets = exclusive_jets(cs; njets = 2, T = EEJet)
                
                # Get jet constituents
                constituent_indices = [constituent_indexes(jet, cs) for jet in jets]
                jet_constituents = ext_mod.build_constituents_cluster(recps, constituent_indices)
                
                # Helper function to compare outputs
                function compare_outputs(actual, reference, name; atol=1e-5, rtol=1e-5)
                    # Handle special float values
                    function normalize_value(x)
                        if isa(x, String)
                            if x == "NaN"
                                return NaN
                            elseif x == "Inf"
                                return Inf
                            elseif x == "-Inf"
                                return -Inf
                            end
                        end
                        return x
                    end
                    
                    # Recursively normalize arrays
                    function normalize_array(arr)
                        if isa(arr, AbstractArray)
                            return [normalize_array(x) for x in arr]
                        else
                            return normalize_value(arr)
                        end
                    end
                    
                    ref_normalized = normalize_array(reference)
                    
                    @testset "$name" begin
                        if isa(actual, Vector{Vector{T}} where T)
                            @test length(actual) == length(ref_normalized)
                            for i in 1:length(actual)
                                @test length(actual[i]) == length(ref_normalized[i])
                                for j in 1:length(actual[i])
                                    act_val = actual[i][j]
                                    ref_val = ref_normalized[i][j]
                                    if isnan(ref_val)
                                        @test isnan(act_val)
                                    elseif isinf(ref_val)
                                        @test isinf(act_val) && sign(act_val) == sign(ref_val)
                                    else
                                        @test isapprox(act_val, ref_val; atol=atol, rtol=rtol)
                                    end
                                end
                            end
                        else
                            @test false  # Unexpected output type
                        end
                    end
                end
                
                @testset "Basic Kinematics" begin
                    compare_outputs(helper_mod.get_pt(jet_constituents), 
                                  reference_outputs["get_pt"], "get_pt")
                    compare_outputs(helper_mod.get_p(jet_constituents), 
                                  reference_outputs["get_p"], "get_p")
                    compare_outputs(helper_mod.get_e(jet_constituents), 
                                  reference_outputs["get_e"], "get_e")
                    compare_outputs(helper_mod.get_type(jet_constituents), 
                                  reference_outputs["get_type"], "get_type")
                    compare_outputs(helper_mod.get_mass(jet_constituents), 
                                  reference_outputs["get_mass"], "get_mass")
                    compare_outputs(helper_mod.get_charge(jet_constituents), 
                                  reference_outputs["get_charge"], "get_charge")
                    compare_outputs(helper_mod.get_theta(jet_constituents), 
                                  reference_outputs["get_theta"], "get_theta")
                    compare_outputs(helper_mod.get_phi(jet_constituents), 
                                  reference_outputs["get_phi"], "get_phi")
                    compare_outputs(helper_mod.get_y(jet_constituents), 
                                  reference_outputs["get_y"], "get_y")
                    compare_outputs(helper_mod.get_eta(jet_constituents), 
                                  reference_outputs["get_eta"], "get_eta")
                    compare_outputs(helper_mod.get_Bz(jet_constituents, tracks), 
                                  reference_outputs["get_Bz"], "get_Bz")
                end
                
                @testset "Track Parameters" begin
                    compare_outputs(helper_mod.get_dxy(jet_constituents, tracks, V, bz), 
                                  reference_outputs["get_dxy"], "get_dxy")
                    compare_outputs(helper_mod.get_dz(jet_constituents, tracks, V, bz), 
                                  reference_outputs["get_dz"], "get_dz")
                    compare_outputs(helper_mod.get_phi0(jet_constituents, tracks, V, bz), 
                                  reference_outputs["get_phi0"], "get_phi0")
                    compare_outputs(helper_mod.get_c(jet_constituents, tracks, bz), 
                                  reference_outputs["get_c"], "get_c")
                    compare_outputs(helper_mod.get_ct(jet_constituents, tracks, bz), 
                                  reference_outputs["get_ct"], "get_ct")
                end
                
                @testset "Covariance Matrix Elements" begin
                    compare_outputs(helper_mod.get_dxydxy(jet_constituents, tracks), 
                                  reference_outputs["get_dxydxy"], "get_dxydxy")
                    compare_outputs(helper_mod.get_dphidxy(jet_constituents, tracks), 
                                  reference_outputs["get_dphidxy"], "get_dphidxy")
                    compare_outputs(helper_mod.get_dphidphi(jet_constituents, tracks), 
                                  reference_outputs["get_dphidphi"], "get_dphidphi")
                    compare_outputs(helper_mod.get_dxyc(jet_constituents, tracks), 
                                  reference_outputs["get_dxyc"], "get_dxyc")
                    compare_outputs(helper_mod.get_phic(jet_constituents, tracks), 
                                  reference_outputs["get_phic"], "get_phic")
                    compare_outputs(helper_mod.get_dptdpt(jet_constituents, tracks), 
                                  reference_outputs["get_dptdpt"], "get_dptdpt")
                    compare_outputs(helper_mod.get_dxydz(jet_constituents, tracks), 
                                  reference_outputs["get_dxydz"], "get_dxydz")
                    compare_outputs(helper_mod.get_phidz(jet_constituents, tracks), 
                                  reference_outputs["get_phidz"], "get_phidz")
                    compare_outputs(helper_mod.get_cdz(jet_constituents, tracks), 
                                  reference_outputs["get_cdz"], "get_cdz")
                    compare_outputs(helper_mod.get_dzdz(jet_constituents, tracks), 
                                  reference_outputs["get_dzdz"], "get_dzdz")
                    compare_outputs(helper_mod.get_dxyctgtheta(jet_constituents, tracks), 
                                  reference_outputs["get_dxyctgtheta"], "get_dxyctgtheta")
                    compare_outputs(helper_mod.get_phictgtheta(jet_constituents, tracks), 
                                  reference_outputs["get_phictgtheta"], "get_phictgtheta")
                    compare_outputs(helper_mod.get_cctgtheta(jet_constituents, tracks), 
                                  reference_outputs["get_cctgtheta"], "get_cctgtheta")
                    compare_outputs(helper_mod.get_dlambdadz(jet_constituents, tracks), 
                                  reference_outputs["get_dlambdadz"], "get_dlambdadz")
                    compare_outputs(helper_mod.get_detadeta(jet_constituents, tracks), 
                                  reference_outputs["get_detadeta"], "get_detadeta")
                end
                
                @testset "Particle Identification" begin
                    compare_outputs(helper_mod.get_is_mu(jet_constituents), 
                                  reference_outputs["get_is_mu"], "get_is_mu")
                    compare_outputs(helper_mod.get_is_el(jet_constituents), 
                                  reference_outputs["get_is_el"], "get_is_el")
                    compare_outputs(helper_mod.get_is_charged_had(jet_constituents), 
                                  reference_outputs["get_is_charged_had"], "get_is_charged_had")
                    compare_outputs(helper_mod.get_is_gamma(jet_constituents), 
                                  reference_outputs["get_is_gamma"], "get_is_gamma")
                    compare_outputs(helper_mod.get_is_neutral_had(jet_constituents), 
                                  reference_outputs["get_is_neutral_had"], "get_is_neutral_had")
                end
                
                @testset "Relative Kinematics" begin
                    compare_outputs(helper_mod.get_erel_cluster(jets, jet_constituents), 
                                  reference_outputs["get_erel_cluster"], "get_erel_cluster")
                    compare_outputs(helper_mod.get_erel_log_cluster(jets, jet_constituents), 
                                  reference_outputs["get_erel_log_cluster"], "get_erel_log_cluster")
                    compare_outputs(helper_mod.get_thetarel_cluster(jets, jet_constituents), 
                                  reference_outputs["get_thetarel_cluster"], "get_thetarel_cluster")
                    compare_outputs(helper_mod.get_phirel_cluster(jets, jet_constituents), 
                                  reference_outputs["get_phirel_cluster"], "get_phirel_cluster")
                    
                    # Combined function returns tuple
                    theta_phi_rel = helper_mod.get_thetarel_phirel_cluster(jets, jet_constituents)
                    compare_outputs(theta_phi_rel[1], 
                                  reference_outputs["get_thetarel_phirel_cluster_theta"], 
                                  "get_thetarel_phirel_cluster_theta")
                    compare_outputs(theta_phi_rel[2], 
                                  reference_outputs["get_thetarel_phirel_cluster_phi"], 
                                  "get_thetarel_phirel_cluster_phi")
                end
                
                @testset "Special Measurements" begin
                    is_charged_had = helper_mod.get_is_charged_had(jet_constituents)
                    compare_outputs(helper_mod.get_dndx(jet_constituents, dNdx, trackdata, is_charged_had), 
                                  reference_outputs["get_dndx"], "get_dndx")
                    compare_outputs(helper_mod.get_mtof(jet_constituents, track_L, trackdata, trackerhits, 
                                                   gammadata, nhdata, calohits, V), 
                                  reference_outputs["get_mtof"], "get_mtof")
                end
                
                @testset "Impact Parameters and Jet Distance" begin
                    # Calculate base values
                    D0 = helper_mod.get_dxy(jet_constituents, tracks, V, bz)
                    Z0 = helper_mod.get_dz(jet_constituents, tracks, V, bz)
                    phi0 = helper_mod.get_phi0(jet_constituents, tracks, V, bz)
                    err2_D0 = helper_mod.get_dxydxy(jet_constituents, tracks)
                    err2_Z0 = helper_mod.get_dzdz(jet_constituents, tracks)
                    
                    # 2D impact parameter
                    sip2d_val = helper_mod.get_Sip2dVal_clusterV(jets, D0, phi0, bz)
                    compare_outputs(sip2d_val, reference_outputs["get_Sip2dVal_clusterV"], 
                                  "get_Sip2dVal_clusterV")
                    
                    btag_sip2d_val = helper_mod.get_btagSip2dVal(jets, D0, phi0, bz)
                    compare_outputs(btag_sip2d_val, reference_outputs["get_btagSip2dVal"], 
                                  "get_btagSip2dVal")
                    
                    compare_outputs(helper_mod.get_Sip2dSig(sip2d_val, err2_D0), 
                                  reference_outputs["get_Sip2dSig"], "get_Sip2dSig")
                    compare_outputs(helper_mod.get_btagSip2dSig(btag_sip2d_val, err2_D0), 
                                  reference_outputs["get_btagSip2dSig"], "get_btagSip2dSig")
                    
                    # 3D impact parameter
                    sip3d_val = helper_mod.get_Sip3dVal_clusterV(jets, D0, Z0, phi0, bz)
                    compare_outputs(sip3d_val, reference_outputs["get_Sip3dVal_clusterV"], 
                                  "get_Sip3dVal_clusterV")
                    
                    btag_sip3d_val = helper_mod.get_btagSip3dVal(jets, D0, Z0, phi0, bz)
                    compare_outputs(btag_sip3d_val, reference_outputs["get_btagSip3dVal"], 
                                  "get_btagSip3dVal")
                    
                    compare_outputs(helper_mod.get_Sip3dSig(sip3d_val, err2_D0, err2_Z0), 
                                  reference_outputs["get_Sip3dSig"], "get_Sip3dSig")
                    compare_outputs(helper_mod.get_btagSip3dSig(btag_sip3d_val, err2_D0, err2_Z0), 
                                  reference_outputs["get_btagSip3dSig"], "get_btagSip3dSig")
                    
                    # Jet distance
                    jet_dist_val = helper_mod.get_JetDistVal_clusterV(jets, jet_constituents, D0, Z0, phi0, bz)
                    compare_outputs(jet_dist_val, reference_outputs["get_JetDistVal_clusterV"], 
                                  "get_JetDistVal_clusterV")
                    
                    btag_jet_dist_val = helper_mod.get_btagJetDistVal(jets, jet_constituents, D0, Z0, phi0, bz)
                    compare_outputs(btag_jet_dist_val, reference_outputs["get_btagJetDistVal"], 
                                  "get_btagJetDistVal")
                    
                    compare_outputs(helper_mod.get_JetDistSig(jet_dist_val, err2_D0, err2_Z0), 
                                  reference_outputs["get_JetDistSig"], "get_JetDistSig")
                    compare_outputs(helper_mod.get_btagJetDistSig(btag_jet_dist_val, err2_D0, err2_Z0), 
                                  reference_outputs["get_btagJetDistSig"], "get_btagJetDistSig")
                end
            end
        end
    end
end