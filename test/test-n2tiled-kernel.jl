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
                              output,
                              event;
                              algorithm = JetAlgorithm.AntiKt,
                              p = nothing,
                              R = 0.4,
                              ptmin = 5.0,
                              recombine = addjets_escheme,
                              preprocess = preprocess_escheme)
    return with_n2tiled_reconstruction(workspace,
                                       event;
                                       algorithm = algorithm,
                                       p = p,
                                       R = R,
                                       recombine = recombine,
                                       preprocess = preprocess) do clusterseq
        inclusive_jets!(output, clusterseq; ptmin = ptmin)
        return nothing
    end
end

function run_owning_n2tiled_event(event)
    return inclusive_jets(tiled_jet_reconstruct(event;
                                                algorithm = JetAlgorithm.AntiKt,
                                                R = 0.4,
                                                preprocess = nothing);
                          ptmin = 5.0)
end

@testset "Reusable N2Tiled reconstruction" begin
    events = read_final_state_particles(events_file_pp)
    event_sizes = length.(events)
    small_event = events[argmin(event_sizes)]
    large_event = events[argmax(event_sizes)]

    @testset "Exact owning/workspace equivalence" begin
        workspace = N2TiledWorkspace()
        output_buffer = LorentzVector{Float64}[]

        cases = ((PseudoJet[], JetAlgorithm.AntiKt, nothing, 0.4),
                 (small_event[1:1], JetAlgorithm.AntiKt, nothing, 0.8),
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

            with_n2tiled_reconstruction(workspace,
                                        event;
                                        algorithm = algorithm,
                                        p = power,
                                        R = distance,
                                        preprocess = nothing) do actual
                @test actual.jets === workspace.jets
                @test actual.history === workspace.history
                test_exact_clustersequence_equality(actual, expected)

                output = inclusive_jets!(output_buffer,
                                         actual;
                                         ptmin = 5.0)
                @test output === output_buffer
                @test output == expected_inclusive
            end
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
    end

    @testset "Owning and borrowed lifetimes" begin
        owning_result = tiled_jet_reconstruct(small_event;
                                              algorithm = JetAlgorithm.AntiKt,
                                              R = 0.4,
                                              preprocess = nothing)
        owning_jets = copy(owning_result.jets)
        owning_history = copy(owning_result.history)

        workspace = N2TiledWorkspace()
        output_buffer = PseudoJet[]
        borrowed_result = Ref{Any}()
        borrowed_output = Ref{Any}()

        with_n2tiled_reconstruction(workspace,
                                    small_event;
                                    algorithm = JetAlgorithm.AntiKt,
                                    R = 0.4,
                                    preprocess = nothing) do clusterseq
            borrowed_result[] = clusterseq
            borrowed_output[] = inclusive_jets!(output_buffer,
                                                clusterseq;
                                                ptmin = 5.0)
        end

        with_n2tiled_reconstruction(workspace,
                                    large_event;
                                    algorithm = JetAlgorithm.AntiKt,
                                    R = 0.4,
                                    preprocess = nothing) do clusterseq
            inclusive_jets!(output_buffer, clusterseq; ptmin = 5.0)
        end

        @test owning_result.jets == owning_jets
        @test owning_result.history == owning_history
        @test borrowed_result[].jets === workspace.jets
        @test borrowed_result[].history === workspace.history
        @test borrowed_output[] === output_buffer
    end

    @testset "In-place inclusive selection" begin
        clusterseq = tiled_jet_reconstruct(small_event;
                                           algorithm = JetAlgorithm.AntiKt,
                                           R = 0.4,
                                           preprocess = nothing)
        output = LorentzVector{Float64}[LorentzVector(1.0, 2.0, 3.0, 4.0)]

        @test inclusive_jets!(output,
                              clusterseq;
                              ptmin = 5.0) === output
        @test output == inclusive_jets(clusterseq; ptmin = 5.0)
        @test_throws ArgumentError inclusive_jets!(clusterseq.jets,
                                                   clusterseq)
    end

    @testset "Capacity reuse and release" begin
        workspace = N2TiledWorkspace()
        output_buffer = LorentzVector{Float64}[]

        run_workspace_event!(workspace, output_buffer, large_event)
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

        run_workspace_event!(workspace, output_buffer, large_event)
        @test workspace.scratch.eta === retained_eta
        @test workspace.scratch.tiledjets === retained_tiledjets
        @test workspace.scratch.NNs === retained_NNs
        @test workspace.scratch.dij === retained_dij
        @test all(workspace.scratch.tiling_cache[shape] === arrays
                  for (shape, arrays) in retained_tilings)

        run_workspace_event!(workspace, output_buffer, small_event)
        @test workspace.jets === retained_jets
        @test workspace.history === retained_history
        @test length(workspace.scratch.eta) == length(small_event)
        @test length(workspace.scratch.tiledjets) == length(small_event)
        @test length(workspace.scratch.NNs) == length(small_event)
        @test length(workspace.scratch.dij) == length(small_event)

        release_n2tiled_workspace_capacity!(workspace)
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

        run_workspace_event!(workspace, output_buffer, small_event)
        @test !isempty(workspace.history)
    end

    @testset "Workspace allocation regression" begin
        event = large_event
        workspace = N2TiledWorkspace()
        output_buffer = LorentzVector{Float64}[]

        run_owning_n2tiled_event(event)
        run_workspace_event!(workspace,
                             output_buffer,
                             event;
                             preprocess = nothing)

        owning_bytes = @allocated run_owning_n2tiled_event(event)
        reused_bytes = @allocated run_workspace_event!(workspace,
                                                       output_buffer,
                                                       event;
                                                       preprocess = nothing)

        @test reused_bytes < owning_bytes ÷ 4
    end
end
