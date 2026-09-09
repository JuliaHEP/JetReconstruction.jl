# Exact owning/workspace output equivalence for every reusable strategy.

include("common.jl")

function test_pp_workspace_output_equivalence!(strategy,
                                               workspace,
                                               event;
                                               algorithm,
                                               p,
                                               R)
    expected_sequence = jet_reconstruct(event;
                                        algorithm = algorithm,
                                        p = p,
                                        R = R,
                                        strategy = strategy,
                                        preprocess = nothing)
    expected_pseudojets = inclusive_jets(expected_sequence,
                                         PseudoJet;
                                         ptmin = 5.0)
    expected_lorentz = inclusive_jets(expected_sequence,
                                      LorentzVector{Float64};
                                      ptmin = 5.0)

    reconstruct = strategy == RecoStrategy.N2Plain ?
                  with_n2plain_reconstruction : with_n2tiled_reconstruction

    reconstruct(workspace,
                event;
                algorithm = algorithm,
                p = p,
                R = R,
                preprocess = nothing) do actual_sequence
        actual_pseudojets = inclusive_jets(actual_sequence,
                                           PseudoJet;
                                           ptmin = 5.0)
        actual_lorentz = inclusive_jets(actual_sequence,
                                        LorentzVector{Float64};
                                        ptmin = 5.0)

        @test actual_pseudojets == expected_pseudojets
        @test actual_lorentz == expected_lorentz
    end

    return nothing
end

@testset "Owning/workspace output equivalence" begin
    @testset "pp N2Plain and N2Tiled" begin
        events = read_final_state_particles(events_file_pp, PseudoJet)
        event_sizes = length.(events)
        event_order = unique(vcat(argmax(event_sizes),
                                  argmin(event_sizes),
                                  collect(eachindex(events))))
        configurations = ((JetAlgorithm.AntiKt, nothing),
                          (JetAlgorithm.CA, nothing),
                          (JetAlgorithm.Kt, nothing),
                          (JetAlgorithm.GenKt, 1.5))

        for strategy in (RecoStrategy.N2Plain, RecoStrategy.N2Tiled)
            workspace = strategy == RecoStrategy.N2Plain ?
                        N2PlainWorkspace(PseudoJet) : N2TiledWorkspace()

            for (algorithm, power) in configurations
                @testset "$(strategy), $(algorithm)" begin
                    for event_index in event_order
                        test_pp_workspace_output_equivalence!(strategy,
                                                              workspace,
                                                              events[event_index];
                                                              algorithm = algorithm,
                                                              p = power,
                                                              R = 0.4)
                    end
                end
            end
        end
    end

    @testset "e+e- N2Plain" begin
        events = read_final_state_particles(events_file_ee, EEJet)
        event_sizes = length.(events)
        small_event = events[argmin(event_sizes)]
        large_event = events[argmax(event_sizes)]
        workspace = N2PlainWorkspace(EEJet)
        cases = ((EEJet[], JetAlgorithm.Durham, nothing, 4.0, nothing),
                 (small_event[1:1], JetAlgorithm.Durham, nothing, 4.0, nothing),
                 (small_event, JetAlgorithm.EEKt, -1, 1.0, nothing),
                 (large_event, JetAlgorithm.EEKt, 1, 2.0, nothing),
                 (large_event, JetAlgorithm.Valencia, 1.2, 0.8, 1.2))

        for (event, algorithm, power, radius, gamma) in cases
            expected_sequence = ee_genkt_algorithm(event;
                                                   algorithm = algorithm,
                                                   p = power,
                                                   R = radius,
                                                   γ = gamma,
                                                   preprocess = nothing)
            expected_eejets = inclusive_jets(expected_sequence,
                                             EEJet;
                                             ptmin = 5.0)
            expected_lorentz = inclusive_jets(expected_sequence,
                                              LorentzVector{Float64};
                                              ptmin = 5.0)

            with_n2plain_reconstruction(workspace,
                                        event;
                                        algorithm = algorithm,
                                        p = power,
                                        R = radius,
                                        γ = gamma,
                                        preprocess = nothing) do actual_sequence
                actual_eejets = inclusive_jets(actual_sequence,
                                               EEJet;
                                               ptmin = 5.0)
                actual_lorentz = inclusive_jets(actual_sequence,
                                                LorentzVector{Float64};
                                                ptmin = 5.0)

                @test actual_eejets == expected_eejets
                @test actual_lorentz == expected_lorentz
            end
        end
    end
end
