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
        inclusive_jets!(workspace, clusterseq; ptmin = ptmin)
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

                output = inclusive_jets!(workspace,
                                         actual;
                                         ptmin = 5.0)
                @test output === workspace.inclusive_output
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

        workspace = N2TiledWorkspace(PseudoJet)
        borrowed_result = Ref{Any}()
        borrowed_output = Ref{Any}()

        with_n2tiled_reconstruction(workspace,
                                    small_event;
                                    algorithm = JetAlgorithm.AntiKt,
                                    R = 0.4,
                                    preprocess = nothing) do clusterseq
            borrowed_result[] = clusterseq
            borrowed_output[] = inclusive_jets!(workspace,
                                                clusterseq;
                                                ptmin = 5.0)
        end

        with_n2tiled_reconstruction(workspace,
                                    large_event;
                                    algorithm = JetAlgorithm.AntiKt,
                                    R = 0.4,
                                    preprocess = nothing) do clusterseq
            inclusive_jets!(workspace, clusterseq; ptmin = 5.0)
        end

        @test owning_result.jets == owning_jets
        @test owning_result.history == owning_history
        @test borrowed_result[].jets === workspace.jets
        @test borrowed_result[].history === workspace.history
        @test borrowed_output[] === workspace.inclusive_output
    end

    @testset "Workspace ownership guard" begin
        workspace = N2TiledWorkspace()

        @test_throws ArgumentError with_n2tiled_reconstruction(workspace,
                                                               small_event;
                                                               algorithm = JetAlgorithm.AntiKt,
                                                               R = 0.4) do _
            with_n2tiled_reconstruction(workspace,
                                        small_event;
                                        algorithm = JetAlgorithm.AntiKt,
                                        R = 0.4) do _
                nothing
            end
        end

        @test_throws ErrorException with_n2tiled_reconstruction(workspace,
                                                                small_event;
                                                                algorithm = JetAlgorithm.AntiKt,
                                                                R = 0.4) do _
            error("intentional callback failure")
        end

        @test with_n2tiled_reconstruction(workspace,
                                          small_event;
                                          algorithm = JetAlgorithm.AntiKt,
                                          R = 0.4) do _
            :recovered
        end == :recovered

        entered = Channel{Nothing}(1)
        release = Channel{Nothing}(1)

        holder = @async with_n2tiled_reconstruction(workspace,
                                                    small_event;
                                                    algorithm = JetAlgorithm.AntiKt,
                                                    R = 0.4) do _
            put!(entered, nothing)
            take!(release)
        end

        wait_status = timedwait(() -> isready(entered) || istaskdone(holder),
                                10.0)
        @test wait_status == :ok

        if isready(entered)
            take!(entered)
            @test_throws ArgumentError run_workspace_event!(workspace,
                                                            small_event)
            put!(release, nothing)
        else
            # Ensure a late holder cannot remain blocked if the wait failed.
            put!(release, nothing)
        end

        fetch(holder)
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

        run_workspace_event!(workspace, large_event)
        retained_jets = workspace.jets
        retained_history = workspace.history
        retained_eta = workspace.scratch.eta
        retained_tiledjets = workspace.scratch.tiledjets
        retained_NNs = workspace.scratch.NNs
        retained_dij = workspace.scratch.dij
        retained_tilings = copy(workspace.scratch.tiling_cache)

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

        release_n2tiled_workspace_capacity!(workspace)
        @test workspace.jets !== retained_jets
        @test workspace.history !== retained_history
        @test workspace.scratch.eta !== retained_eta
        @test workspace.scratch.tiledjets !== retained_tiledjets
        @test workspace.scratch.NNs !== retained_NNs
        @test workspace.scratch.dij !== retained_dij
        @test isempty(workspace.jets)
        @test isempty(workspace.history)
        @test isempty(workspace.inclusive_output)
        @test isempty(workspace.scratch.eta)
        @test isempty(workspace.scratch.tiledjets)
        @test isempty(workspace.scratch.NNs)
        @test isempty(workspace.scratch.dij)
        @test isempty(workspace.scratch.tiling_cache)

        run_workspace_event!(workspace, small_event)
        @test !isempty(workspace.history)
    end

    @testset "Workspace allocation regression" begin
        event = large_event
        workspace = N2TiledWorkspace()

        run_owning_n2tiled_event(event)
        run_workspace_event!(workspace, event; preprocess = nothing)

        owning_bytes = @allocated run_owning_n2tiled_event(event)
        reused_bytes = @allocated run_workspace_event!(workspace,
                                                       event;
                                                       preprocess = nothing)

        @test reused_bytes < owning_bytes ÷ 4
    end

    @testset "Abstract-vector tiling input" begin
        eta_values = [-2.0, -0.5, 0.25, 1.5]
        eta_view = @view eta_values[1:3]

        @test JetReconstruction.determine_rapidity_extent(eta_view) ==
              JetReconstruction.determine_rapidity_extent(collect(eta_view))
        @test JetReconstruction.setup_tiling(eta_view, 0.4) ==
              JetReconstruction.setup_tiling(collect(eta_view), 0.4)
    end
end
