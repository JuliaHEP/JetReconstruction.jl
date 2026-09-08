# Regression tests for reusable N2Tiled reconstruction storage.

include("common.jl")

function test_exact_clustersequence_equality(actual, expected)
    @test actual.algorithm == expected.algorithm
    @test actual.power == expected.power
    @test actual.R == expected.R
    @test actual.strategy == expected.strategy
    @test actual.n_initial_jets == expected.n_initial_jets
    @test actual.Qtot == expected.Qtot
    @test actual.jets == expected.jets
    @test actual.history == expected.history

    return nothing
end

function run_workspace_event!(workspace,
                              event;
                              algorithm = JetAlgorithm.AntiKt,
                              p = nothing,
                              R = 0.4,
                              recombine = addjets_escheme,
                              preprocess = preprocess_escheme)
    return with_n2tiled_reconstruction(workspace,
                                       event;
                                       algorithm = algorithm,
                                       p = p,
                                       R = R,
                                       recombine = recombine,
                                       preprocess = preprocess) do _
        return nothing
    end
end

function run_owning_n2tiled_event(event)
    return tiled_jet_reconstruct(event;
                                 algorithm = JetAlgorithm.AntiKt,
                                 R = 0.4,
                                 preprocess = nothing)
end

@testset "Reusable N2Tiled reconstruction" begin
    @test :release_n2tiled_workspace_capacity! ∉ names(JetReconstruction)

    events = read_final_state_particles(events_file_pp)
    event_sizes = length.(events)
    small_event = events[argmin(event_sizes)]
    large_event = events[argmax(event_sizes)]

    @testset "Exact owning/workspace equivalence" begin
        workspace = N2TiledWorkspace()

        cases = ((small_event[1:1], JetAlgorithm.AntiKt, nothing, 0.8),
                 (small_event, JetAlgorithm.AntiKt, nothing, 0.4),
                 (large_event, JetAlgorithm.AntiKt, nothing, 0.2),
                 (small_event, JetAlgorithm.Kt, nothing, 0.4),
                 (large_event, JetAlgorithm.CA, nothing, 0.8),
                 (small_event, JetAlgorithm.GenKt, 1.5, 0.4))

        for (event, algorithm, power, distance) in cases
            expected = tiled_jet_reconstruct(event;
                                             algorithm = algorithm,
                                             p = power,
                                             R = distance,
                                             preprocess = nothing)
            expected_inclusive = inclusive_jets(expected; ptmin = 5.0)

            output = with_n2tiled_reconstruction(workspace,
                                                 event;
                                                 algorithm = algorithm,
                                                 p = power,
                                                 R = distance,
                                                 preprocess = nothing) do actual
                @test actual.jets === workspace.jets
                @test actual.history === workspace.history
                @test length(actual.history) == 2 * length(event)
                test_exact_clustersequence_equality(actual, expected)

                return inclusive_jets(actual; ptmin = 5.0)
            end
            @test output == expected_inclusive
        end
    end

    @testset "Preprocessing and recombination compatibility" begin
        workspace = N2TiledWorkspace()
        event = small_event

        configurations = ((preprocess_escheme, addjets_escheme),
                          (preprocess_ptscheme, addjets_ptscheme),
                          (preprocess_pt2scheme, addjets_pt2scheme))

        for (preprocess, recombine) in configurations
            expected = tiled_jet_reconstruct(event;
                                             algorithm = JetAlgorithm.AntiKt,
                                             R = 0.4,
                                             preprocess = preprocess,
                                             recombine = recombine)

            with_n2tiled_reconstruction(workspace,
                                        event;
                                        algorithm = JetAlgorithm.AntiKt,
                                        R = 0.4,
                                        preprocess = preprocess,
                                        recombine = recombine) do actual
                test_exact_clustersequence_equality(actual, expected)
            end
        end

        lorentz_event = lorentzvector.(event)
        expected = tiled_jet_reconstruct(lorentz_event;
                                         algorithm = JetAlgorithm.AntiKt,
                                         R = 0.4)

        with_n2tiled_reconstruction(workspace,
                                    lorentz_event;
                                    algorithm = JetAlgorithm.AntiKt,
                                    R = 0.4) do actual
            test_exact_clustersequence_equality(actual, expected)
        end

        expected = tiled_jet_reconstruct(lorentz_event;
                                         algorithm = JetAlgorithm.AntiKt,
                                         R = 0.4,
                                         preprocess = nothing)

        with_n2tiled_reconstruction(workspace,
                                    lorentz_event;
                                    algorithm = JetAlgorithm.AntiKt,
                                    R = 0.4,
                                    preprocess = nothing) do actual
            test_exact_clustersequence_equality(actual, expected)
        end
    end

    @testset "Owning and borrowed lifetimes" begin
        owning_result = tiled_jet_reconstruct(small_event;
                                              algorithm = JetAlgorithm.AntiKt,
                                              R = 0.4,
                                              preprocess = nothing)
        owning_jets = copy(owning_result.jets)
        owning_history = copy(owning_result.history)

        workspace = N2TiledWorkspace()
        borrowed_result = Ref{Any}()
        owned_output = Ref{Any}()
        retained_result = Ref{Any}()

        with_n2tiled_reconstruction(workspace,
                                    small_event;
                                    algorithm = JetAlgorithm.AntiKt,
                                    R = 0.4,
                                    preprocess = nothing) do clusterseq
            borrowed_result[] = clusterseq
            owned_output[] = inclusive_jets(clusterseq,
                                            PseudoJet;
                                            ptmin = 5.0)
            retained_result[] = deepcopy(clusterseq)
        end

        owned_output_snapshot = copy(owned_output[])

        with_n2tiled_reconstruction(workspace,
                                    large_event;
                                    algorithm = JetAlgorithm.AntiKt,
                                    R = 0.4,
                                    preprocess = nothing) do _
            nothing
        end

        @test owning_result.jets == owning_jets
        @test owning_result.history == owning_history
        @test borrowed_result[].jets === workspace.jets
        @test borrowed_result[].history === workspace.history
        @test owned_output[] == owned_output_snapshot
        @test retained_result[].jets == owning_jets
        @test retained_result[].history == owning_history
        @test retained_result[].jets !== workspace.jets
        @test retained_result[].history !== workspace.history
    end

    @testset "Capacity reuse and release" begin
        workspace = N2TiledWorkspace()

        run_workspace_event!(workspace, large_event)
        retained_jets = workspace.jets
        retained_history = workspace.history
        retained_eta = workspace.scratch.eta
        retained_tiledjets = workspace.scratch.tiledjets
        retained_NNs = workspace.scratch.NNs
        retained_dij = workspace.scratch.dij
        retained_tilings = copy(workspace.scratch.tiling_cache)

        @test length(workspace.scratch.eta) == length(large_event)
        @test length(workspace.scratch.tiledjets) == length(large_event)
        @test length(workspace.scratch.NNs) == length(large_event)
        @test length(workspace.scratch.dij) == length(large_event)

        run_workspace_event!(workspace, large_event)
        @test workspace.scratch.eta === retained_eta
        @test workspace.scratch.tiledjets === retained_tiledjets
        @test workspace.scratch.NNs === retained_NNs
        @test workspace.scratch.dij === retained_dij
        @test all(workspace.scratch.tiling_cache[shape] === arrays
                  for (shape, arrays) in retained_tilings)

        run_workspace_event!(workspace, small_event)
        @test workspace.jets === retained_jets
        @test workspace.history === retained_history
        @test length(workspace.scratch.eta) == length(small_event)
        @test length(workspace.scratch.tiledjets) == length(small_event)
        @test length(workspace.scratch.NNs) == length(small_event)
        @test length(workspace.scratch.dij) == length(small_event)

        JetReconstruction.release_n2tiled_workspace_capacity!(workspace)
        @test workspace.jets !== retained_jets
        @test workspace.history !== retained_history
        @test workspace.scratch.eta !== retained_eta
        @test workspace.scratch.tiledjets !== retained_tiledjets
        @test workspace.scratch.NNs !== retained_NNs
        @test workspace.scratch.dij !== retained_dij
        @test isempty(workspace.jets)
        @test isempty(workspace.history)
        @test isempty(workspace.scratch.eta)
        @test isempty(workspace.scratch.tiledjets)
        @test isempty(workspace.scratch.NNs)
        @test isempty(workspace.scratch.dij)
        @test isempty(workspace.scratch.tiling_cache)

        run_workspace_event!(workspace, small_event)
        @test !isempty(workspace.history)
    end

    @testset "Workspace allocation regression" begin
        @test_throws ErrorException JetReconstruction._sizehint_for_reuse!(Int[], -1)

        event = large_event
        workspace = N2TiledWorkspace()

        run_owning_n2tiled_event(event)
        run_workspace_event!(workspace,
                             event;
                             preprocess = nothing)

        owning_bytes = @allocated run_owning_n2tiled_event(event)
        reused_bytes = @allocated run_workspace_event!(workspace,
                                                       event;
                                                       preprocess = nothing)

        @test reused_bytes < owning_bytes ÷ 4
    end
end
