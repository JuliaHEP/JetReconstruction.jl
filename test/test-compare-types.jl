# Compare inputting data in PseudoJet with using a LorentzVector

include("common.jl")

function do_test_compare_types(strategy::RecoStrategy.Strategy;
                               algorithm::JetAlgorithm.Algorithm,
                               ptmin::Float64 = 5.0,
                               distance::Float64 = 0.4,
                               power::Union{Real, Nothing} = nothing)

    # Strategy
    if (strategy == RecoStrategy.N2Plain)
        jet_reconstruction = plain_jet_reconstruct
        strategy_name = "N2Plain"
    elseif (strategy == RecoStrategy.N2Tiled)
        jet_reconstruction = tiled_jet_reconstruct
        strategy_name = "N2Tiled"
    else
        throw(ErrorException("Strategy not yet implemented"))
    end

    # Now run our jet reconstruction...
    # From PseudoJets
    events::Vector{Vector{PseudoJet}} = read_final_state_particles(events_file_pp)
    jet_collection = FinalJets[]
    for (ievt, event) in enumerate(events)
        finaljets = final_jets(inclusive_jets(jet_reconstruction(event, R = distance,
                                                                 algorithm = algorithm,
                                                                 p = power), ptmin = ptmin))
        sort_jets!(finaljets)
        push!(jet_collection, FinalJets(ievt, finaljets))
    end

    # From LorentzVector
    events_lv = read_final_state_particles(events_file_pp, LorentzVector{Float64})
    jet_collection_lv = FinalJets[]
    for (ievt, event) in enumerate(events_lv)
        finaljets = final_jets(inclusive_jets(jet_reconstruction(event, R = distance,
                                                                 algorithm = algorithm,
                                                                 p = power), ptmin = ptmin))
        sort_jets!(finaljets)
        push!(jet_collection_lv, FinalJets(ievt, finaljets))
    end

    @testset "Jet Reconstruction Compare PseudoJet and LorentzVector, Strategy $strategy_name, Algorithm $algorithm" begin
        # Here we test that inputting LorentzVector gave the same results as PseudoJets
        for (ievt, (event, event_lv)) in enumerate(zip(jet_collection, jet_collection_lv))
            @testset "Event $(ievt)" begin
                @test size(event.jets) == size(event_lv.jets)
                # Test each jet in turn
                for (jet, jet_lv) in zip(event.jets, event_lv.jets)
                    @test jet == jet_lv
                end
            end
        end
    end
end

do_test_compare_types(RecoStrategy.N2Plain, algorithm = JetAlgorithm.AntiKt)
do_test_compare_types(RecoStrategy.N2Tiled, algorithm = JetAlgorithm.AntiKt)
