# Regression tests for reusable N2Plain reconstruction storage.

include("common.jl")

function test_n2plain_clustersequence_equality(actual, expected)
    @test actual.algorithm == expected.algorithm
    @test actual.power == expected.power
    @test actual.R == expected.R
    @test actual.strategy == RecoStrategy.N2Plain
    @test actual.strategy == expected.strategy
    @test actual.n_initial_jets == expected.n_initial_jets
    @test actual.Qtot == expected.Qtot
    @test actual.jets == expected.jets
    @test actual.history == expected.history

    return nothing
end

function run_pp_n2plain_workspace!(workspace,
                                   event;
                                   algorithm = JetAlgorithm.AntiKt,
                                   p = nothing,
                                   R = 0.4,
                                   preprocess = nothing,
                                   recombine = addjets_escheme)
    return with_n2plain_reconstruction(workspace,
                                       event;
                                       algorithm = algorithm,
                                       p = p,
                                       R = R,
                                       preprocess = preprocess,
                                       recombine = recombine) do clusterseq
        return inclusive_jets(clusterseq; ptmin = 5.0)
    end
end

function run_ee_n2plain_workspace!(workspace,
                                   event;
                                   algorithm = JetAlgorithm.Durham,
                                   p = nothing,
                                   R = 4.0,
                                   γ = nothing,
                                   preprocess = nothing,
                                   recombine = addjets_escheme)
    return with_n2plain_reconstruction(workspace,
                                       event;
                                       algorithm = algorithm,
                                       p = p,
                                       R = R,
                                       γ = γ,
                                       preprocess = preprocess,
                                       recombine = recombine) do clusterseq
        return inclusive_jets(clusterseq; ptmin = 5.0)
    end
end

function run_pp_n2plain_reconstruction_only!(workspace, event)
    return with_n2plain_reconstruction(workspace,
                                       event;
                                       algorithm = JetAlgorithm.AntiKt,
                                       R = 0.4,
                                       preprocess = nothing) do _
        return nothing
    end
end

function pp_scratch_vectors(workspace)
    scratch = workspace.scratch

    return (scratch.kt2,
            scratch.phi,
            scratch.rapidity,
            scratch.nn,
            scratch.nndist,
            scratch.nndij,
            scratch.clusterseq_index)
end

function ee_scratch_vectors(workspace)
    scratch = workspace.scratch

    return (scratch.index,
            scratch.nni,
            scratch.nndist,
            scratch.dijdist,
            scratch.nx,
            scratch.ny,
            scratch.nz,
            scratch.E2p)
end

@testset "Reusable N2Plain reconstruction" begin
    @test :release_n2plain_workspace_capacity! ∉ names(JetReconstruction)

    pp_events = read_final_state_particles(events_file_pp, PseudoJet)
    ee_events = read_final_state_particles(events_file_ee, EEJet)
    pp_sizes = length.(pp_events)
    ee_sizes = length.(ee_events)
    small_pp = pp_events[argmin(pp_sizes)]
    large_pp = pp_events[argmax(pp_sizes)]
    small_ee = ee_events[argmin(ee_sizes)]
    large_ee = ee_events[argmax(ee_sizes)]

    @testset "pp owning/workspace ClusterSequence equivalence" begin
        workspace = N2PlainWorkspace(PseudoJet)
        cases = ((PseudoJet[], JetAlgorithm.AntiKt, nothing, 0.4),
                 (small_pp[1:1], JetAlgorithm.Kt, nothing, 0.8),
                 (small_pp, JetAlgorithm.AntiKt, nothing, 0.4),
                 (large_pp, JetAlgorithm.CA, nothing, 0.8),
                 (small_pp, JetAlgorithm.GenKt, -1.0, 0.4),
                 (large_pp, JetAlgorithm.GenKt, 1.5, 0.4))

        for (event, algorithm, power, radius) in cases
            expected = plain_jet_reconstruct(event;
                                             algorithm = algorithm,
                                             p = power,
                                             R = radius,
                                             preprocess = nothing)
            expected_inclusive = inclusive_jets(expected; ptmin = 5.0)

            with_n2plain_reconstruction(workspace,
                                        event;
                                        algorithm = algorithm,
                                        p = power,
                                        R = radius,
                                        preprocess = nothing) do actual
                @test actual.jets === workspace.jets
                @test actual.history === workspace.history
                test_n2plain_clustersequence_equality(actual, expected)

                @test inclusive_jets(actual; ptmin = 5.0) == expected_inclusive
            end
        end
    end

    @testset "e+e- owning/workspace ClusterSequence equivalence" begin
        workspace = N2PlainWorkspace(EEJet)
        cases = ((EEJet[], JetAlgorithm.Durham, nothing, 4.0, nothing),
                 (small_ee[1:1], JetAlgorithm.Durham, nothing, 4.0, nothing),
                 (small_ee, JetAlgorithm.Durham, nothing, 4.0, nothing),
                 (large_ee, JetAlgorithm.EEKt, -1, 1.0, nothing),
                 (small_ee, JetAlgorithm.EEKt, 1, 2.0, nothing),
                 (large_ee, JetAlgorithm.Valencia, 1.2, 0.8, 1.2))

        for (event, algorithm, power, radius, gamma) in cases
            expected = ee_genkt_algorithm(event;
                                          algorithm = algorithm,
                                          p = power,
                                          R = radius,
                                          γ = gamma,
                                          preprocess = nothing)
            expected_inclusive = inclusive_jets(expected; ptmin = 5.0)

            with_n2plain_reconstruction(workspace,
                                        event;
                                        algorithm = algorithm,
                                        p = power,
                                        R = radius,
                                        γ = gamma,
                                        preprocess = nothing) do actual
                @test actual.jets === workspace.jets
                @test actual.history === workspace.history
                test_n2plain_clustersequence_equality(actual, expected)

                @test inclusive_jets(actual; ptmin = 5.0) == expected_inclusive
            end
        end
    end

    @testset "Preprocessing and recombination compatibility" begin
        workspace = N2PlainWorkspace(PseudoJet)
        configurations = ((preprocess_escheme, addjets_escheme),
                          (preprocess_ptscheme, addjets_ptscheme),
                          (preprocess_pt2scheme, addjets_pt2scheme))

        for (preprocess, recombine) in configurations
            expected = plain_jet_reconstruct(small_pp;
                                             algorithm = JetAlgorithm.AntiKt,
                                             R = 0.4,
                                             preprocess = preprocess,
                                             recombine = recombine)

            with_n2plain_reconstruction(workspace,
                                        small_pp;
                                        algorithm = JetAlgorithm.AntiKt,
                                        R = 0.4,
                                        preprocess = preprocess,
                                        recombine = recombine) do actual
                test_n2plain_clustersequence_equality(actual, expected)
            end
        end

        lorentz_event = lorentzvector.(small_ee)
        ee_workspace = N2PlainWorkspace(EEJet)
        expected = ee_genkt_algorithm(lorentz_event;
                                      algorithm = JetAlgorithm.Durham)

        with_n2plain_reconstruction(ee_workspace,
                                    lorentz_event;
                                    algorithm = JetAlgorithm.Durham) do actual
            test_n2plain_clustersequence_equality(actual, expected)
        end
    end

    @testset "Borrowed full ClusterSequence queries" begin
        event = first(pp_events)
        workspace = N2PlainWorkspace(PseudoJet)
        expected = plain_jet_reconstruct(event;
                                         algorithm = JetAlgorithm.Kt,
                                         R = 1.0,
                                         preprocess = nothing)
        expected_exclusive = exclusive_jets(expected, PseudoJet; njets = 4)
        expected_inclusive = inclusive_jets(expected, PseudoJet; ptmin = 5.0)

        with_n2plain_reconstruction(workspace,
                                    event;
                                    algorithm = JetAlgorithm.Kt,
                                    R = 1.0,
                                    preprocess = nothing) do actual
            @test exclusive_jets(actual, PseudoJet; njets = 4) == expected_exclusive

            actual_inclusive = inclusive_jets(actual, PseudoJet; ptmin = 5.0)
            @test actual_inclusive == expected_inclusive
            @test constituent_indexes(first(actual_inclusive), actual) ==
                  constituent_indexes(first(expected_inclusive), expected)
            @test parent_jets(first(actual_inclusive), actual) ==
                  parent_jets(first(expected_inclusive), expected)
        end
    end

    @testset "Workspace family validation" begin
        pp_workspace = N2PlainWorkspace(PseudoJet)
        ee_workspace = N2PlainWorkspace(EEJet)

        @test_throws ArgumentError run_pp_n2plain_workspace!(ee_workspace, small_pp)
        @test_throws ArgumentError run_ee_n2plain_workspace!(pp_workspace, small_ee)
    end

    @testset "Exact scratch lengths, reuse, and release" begin
        pp_workspace = N2PlainWorkspace(PseudoJet)

        run_pp_n2plain_workspace!(pp_workspace, large_pp)
        retained_pp_jets = pp_workspace.jets
        retained_pp_history = pp_workspace.history
        retained_pp_scratch = pp_scratch_vectors(pp_workspace)
        @test all(length(array) == length(large_pp) for array in retained_pp_scratch)

        run_pp_n2plain_workspace!(pp_workspace, small_pp)
        @test pp_workspace.jets === retained_pp_jets
        @test pp_workspace.history === retained_pp_history
        @test all(current === retained
                  for (current, retained) in zip(pp_scratch_vectors(pp_workspace),
                                                 retained_pp_scratch))
        @test all(length(array) == length(small_pp)
                  for array in pp_scratch_vectors(pp_workspace))

        JetReconstruction.release_n2plain_workspace_capacity!(pp_workspace)
        @test pp_workspace.jets !== retained_pp_jets
        @test pp_workspace.history !== retained_pp_history
        @test all(current !== retained
                  for (current, retained) in zip(pp_scratch_vectors(pp_workspace),
                                                 retained_pp_scratch))
        @test isempty(pp_workspace.jets)
        @test isempty(pp_workspace.history)
        @test all(isempty, pp_scratch_vectors(pp_workspace))

        ee_workspace = N2PlainWorkspace(EEJet)

        run_ee_n2plain_workspace!(ee_workspace, large_ee)
        retained_ee_jets = ee_workspace.jets
        retained_ee_history = ee_workspace.history
        retained_ee_scratch_object = ee_workspace.scratch
        retained_ee_scratch = ee_scratch_vectors(ee_workspace)
        @test all(length(array) == length(large_ee) for array in retained_ee_scratch)

        run_ee_n2plain_workspace!(ee_workspace, small_ee)
        @test ee_workspace.jets === retained_ee_jets
        @test ee_workspace.history === retained_ee_history
        @test ee_workspace.scratch === retained_ee_scratch_object
        @test all(current === retained
                  for (current, retained) in zip(ee_scratch_vectors(ee_workspace),
                                                 retained_ee_scratch))
        @test all(length(array) == length(small_ee)
                  for array in ee_scratch_vectors(ee_workspace))

        JetReconstruction.release_n2plain_workspace_capacity!(ee_workspace)
        @test ee_workspace.jets !== retained_ee_jets
        @test ee_workspace.history !== retained_ee_history
        @test ee_workspace.scratch !== retained_ee_scratch_object
        @test isempty(ee_workspace.jets)
        @test isempty(ee_workspace.history)
        @test all(isempty, ee_scratch_vectors(ee_workspace))
    end

    @testset "Owning and borrowed lifetimes" begin
        owning_result = plain_jet_reconstruct(small_pp;
                                              algorithm = JetAlgorithm.AntiKt,
                                              R = 0.4,
                                              preprocess = nothing)
        owning_jets = copy(owning_result.jets)
        owning_history = copy(owning_result.history)
        workspace = N2PlainWorkspace(PseudoJet)
        borrowed_result = Ref{Any}()
        borrowed_output = Ref{Any}()

        with_n2plain_reconstruction(workspace,
                                    small_pp;
                                    algorithm = JetAlgorithm.AntiKt,
                                    R = 0.4,
                                    preprocess = nothing) do clusterseq
            borrowed_result[] = clusterseq
            borrowed_output[] = inclusive_jets(clusterseq; ptmin = 5.0)
        end

        run_pp_n2plain_workspace!(workspace, large_pp)

        @test owning_result.jets == owning_jets
        @test owning_result.history == owning_history
        @test borrowed_result[].jets === workspace.jets
        @test borrowed_result[].history === workspace.history
        @test borrowed_output[] isa Vector{LorentzVector{Float64}}
    end

    @testset "Workspace allocation regression" begin
        pp_workspace = N2PlainWorkspace(PseudoJet)

        run_pp_n2plain_reconstruction_only!(pp_workspace, large_pp)
        plain_jet_reconstruct(large_pp;
                              algorithm = JetAlgorithm.AntiKt,
                              R = 0.4,
                              preprocess = nothing)

        reused_bytes = @allocated run_pp_n2plain_reconstruction_only!(pp_workspace,
                                                                      large_pp)
        owning_bytes = @allocated plain_jet_reconstruct(large_pp;
                                                        algorithm = JetAlgorithm.AntiKt,
                                                        R = 0.4,
                                                        preprocess = nothing)
        @test reused_bytes < owning_bytes ÷ 4
    end

    @testset "One workspace per persistent task" begin
        sample_events = pp_events[1:min(12, length(pp_events))]
        expected = [inclusive_jets(plain_jet_reconstruct(event;
                                                         algorithm = JetAlgorithm.AntiKt,
                                                         R = 0.4,
                                                         preprocess = nothing);
                                   ptmin = 5.0) for event in sample_events]
        results = Vector{Vector{LorentzVector{Float64}}}(undef, length(sample_events))
        next_index = Threads.Atomic{Int}(1)
        worker_count = min(max(Threads.nthreads(), 1), 4)

        @sync for _ in 1:worker_count
            Threads.@spawn begin
                workspace = N2PlainWorkspace(PseudoJet)

                while true
                    event_index = Threads.atomic_add!(next_index, 1)
                    event_index > length(sample_events) && break

                    with_n2plain_reconstruction(workspace,
                                                sample_events[event_index];
                                                algorithm = JetAlgorithm.AntiKt,
                                                R = 0.4,
                                                preprocess = nothing) do clusterseq
                        results[event_index] = inclusive_jets(clusterseq;
                                                              ptmin = 5.0)
                    end
                end
            end
        end

        @test results == expected
    end
end
