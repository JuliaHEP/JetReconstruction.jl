# Tests functions exposed as C-interface

include("common.jl")
using JetReconstruction.C_JetReconstruction

# helper function make generic field-wise comparison
function struct_approx_equal(x::T, y::T) where {T}
    return all(getfield(x, f) ≈ getfield(y, f) for f in fieldnames(typeof(x)))
end

function compare_results(ptr::Ptr{C_JetReconstruction.C_ClusterSequence{T}},
                         cluster_seq::ClusterSequence{T}) where {T}
    @test ptr != C_NULL
    c_cluster_seq = unsafe_load(ptr)
    @test c_cluster_seq.algorithm == cluster_seq.algorithm
    @test c_cluster_seq.power ≈ cluster_seq.power
    @test c_cluster_seq.R ≈ cluster_seq.R
    @test c_cluster_seq.strategy == cluster_seq.strategy
    @test c_cluster_seq.jets_length == length(cluster_seq.jets)
    c_jets = C_JetReconstruction.unsafe_wrap_c_array(c_cluster_seq.jets,
                                                     c_cluster_seq.jets_length)
    @test all(struct_approx_equal.(c_jets, cluster_seq.jets))
    @test c_cluster_seq.n_initial_jets == cluster_seq.n_initial_jets
    @test c_cluster_seq.history_length == length(cluster_seq.history)
    c_history = C_JetReconstruction.unsafe_wrap_c_array(c_cluster_seq.history,
                                                        c_cluster_seq.history_length)
    @test all(struct_approx_equal.(c_history, cluster_seq.history))
    @test c_cluster_seq.Qtot ≈ cluster_seq.Qtot
end

function compare_results(ptr::Ptr{C_JetReconstruction.C_JetsResult{T}},
                         jets::Vector{T}) where {T}
    @test ptr != C_NULL
    c_results = unsafe_load(ptr)
    @test c_results.length == length(jets)
    c_data = C_JetReconstruction.unsafe_wrap_c_array(c_results.data,
                                                     c_results.length)
    @test all(struct_approx_equal.(c_data, jets))
end

function test_pseudojet()
    ptr = Ptr{PseudoJet}(Libc.malloc(sizeof(PseudoJet)))
    @test ptr != C_NULL
    ret = C_JetReconstruction.jetreconstruction_PseudoJet_init(ptr, 0.1, 0.2, 0.3, 1.0, 1)
    @test C_JetReconstruction.StatusCode.T(ret) == C_JetReconstruction.StatusCode.OK
    c_jet = unsafe_load(ptr)
    jet = PseudoJet(0.1, 0.2, 0.3, 1.0; cluster_hist_index = 1)
    struct_approx_equal(jet, c_jet)
end

function test_jet_reconstruct(filename; algorithm, R, strategy, power = nothing,
                              T = PseudoJet)
    @testset "C-interface jet reconstruct" begin
        events = JetReconstruction.read_final_state_particles(filename)
        results = Vector{ClusterSequence{T}}(undef, length(events))
        c_results = Vector{Ptr{C_JetReconstruction.C_ClusterSequence{T}}}(undef,
                                                                          length(events))
        for (ievent, event) in enumerate(events)
            c_event, c_event_length = C_JetReconstruction.make_c_array(event)
            cluster_seq_ptr = Ptr{C_JetReconstruction.C_ClusterSequence{T}}(Libc.malloc(sizeof(C_JetReconstruction.C_ClusterSequence{T})))
            @test cluster_seq_ptr != C_NULL

            ret = C_JetReconstruction.jetreconstruction_jet_reconstruct(c_event,
                                                                        c_event_length,
                                                                        algorithm,
                                                                        R,
                                                                        strategy,
                                                                        cluster_seq_ptr)
            @test C_JetReconstruction.StatusCode.T(ret) == C_JetReconstruction.StatusCode.OK

            cluster_seq = JetReconstruction.jet_reconstruct(event; R = R,
                                                            p = power,
                                                            algorithm = algorithm,
                                                            strategy = strategy)
            compare_results(cluster_seq_ptr, cluster_seq)
            @inbounds results[ievent] = cluster_seq
            @inbounds c_results[ievent] = cluster_seq_ptr
            Libc.free(c_event)
        end
        return c_results, results
    end
end

function test_inclusive_jets(cluster_seq_ptrs::Vector{Ptr{C_JetReconstruction.C_ClusterSequence{T}}},
                             cluster_seqs::Vector{JetReconstruction.ClusterSequence{T}};
                             ptmin) where {T}
    @testset "C-interface inclusive jets" begin
        results_ptr = Ptr{C_JetReconstruction.C_JetsResult{T}}(Libc.malloc(sizeof(C_JetReconstruction.C_JetsResult{T})))
        for (cluster_seq_ptr, cluster_seq) in zip(cluster_seq_ptrs, cluster_seqs)
            ret = C_JetReconstruction.jetreconstruction_inclusive_jets(cluster_seq_ptr,
                                                                       ptmin,
                                                                       results_ptr)
            @test C_JetReconstruction.StatusCode.T(ret) == C_JetReconstruction.StatusCode.OK
            results = inclusive_jets(cluster_seq; ptmin = ptmin, T = T)
            compare_results(results_ptr, results)
            C_JetReconstruction.jetreconstruction_JetsResult_free_members_(results_ptr)
        end
        Libc.free(results_ptr)
    end
end

function test_exclusive_jets_njets(cluster_seq_ptrs::Vector{Ptr{C_JetReconstruction.C_ClusterSequence{T}}},
                                   cluster_seqs::Vector{JetReconstruction.ClusterSequence{T}};
                                   njets) where {T}
    @testset "C-interface exclusive jets njets" begin
        results_ptr = Ptr{C_JetReconstruction.C_JetsResult{T}}(Libc.malloc(sizeof(C_JetReconstruction.C_JetsResult{T})))
        for (cluster_seq_ptr, cluster_seq) in zip(cluster_seq_ptrs, cluster_seqs)
            ret = C_JetReconstruction.jetreconstruction_exclusive_jets_njets(cluster_seq_ptr,
                                                                             Csize_t(njets),
                                                                             results_ptr)
            @test C_JetReconstruction.StatusCode.T(ret) == C_JetReconstruction.StatusCode.OK
            results = exclusive_jets(cluster_seq; njets = njets, T = T)
            compare_results(results_ptr, results)
            C_JetReconstruction.jetreconstruction_JetsResult_free_members_(results_ptr)
        end
        Libc.free(results_ptr)
    end
end

function test_exclusive_jets_dcut(cluster_seq_ptrs::Vector{Ptr{C_JetReconstruction.C_ClusterSequence{T}}},
                                  cluster_seqs::Vector{JetReconstruction.ClusterSequence{T}};
                                  dcut) where {T}
    @testset "C-interface exclusive jets dcut" begin
        results_ptr = Ptr{C_JetReconstruction.C_JetsResult{T}}(Libc.malloc(sizeof(C_JetReconstruction.C_JetsResult{T})))
        for (cluster_seq_ptr, cluster_seq) in zip(cluster_seq_ptrs, cluster_seqs)
            ret = C_JetReconstruction.jetreconstruction_exclusive_jets_dcut(cluster_seq_ptr,
                                                                            dcut,
                                                                            results_ptr)
            @test C_JetReconstruction.StatusCode.T(ret) == C_JetReconstruction.StatusCode.OK
            results = exclusive_jets(cluster_seq; dcut = dcut, T = T)
            compare_results(results_ptr, results)
            C_JetReconstruction.jetreconstruction_JetsResult_free_members_(results_ptr)
        end
        Libc.free(results_ptr)
    end
end

@testset "C-interface JetReconstruction pp" begin
    test_cone_size = 0.4

    test_pseudojet()

    for alg in [JetAlgorithm.AntiKt, JetAlgorithm.CA, JetAlgorithm.Kt],
        stg in [RecoStrategy.Best, RecoStrategy.N2Plain, RecoStrategy.N2Tiled]

        power = JetReconstruction.algorithm2power[alg]
        @testset "C-interface JetReconstruction comparison: alg=$alg, p=$power, R=$test_cone_size, strategy=$stg" begin
            cluster_seq_ptrs, cluster_seqs = test_jet_reconstruct(events_file_pp;
                                                                  algorithm = alg,
                                                                  strategy = stg,
                                                                  R = test_cone_size,
                                                                  power = power)
            test_inclusive_jets(cluster_seq_ptrs, cluster_seqs; ptmin = 5.0)
            if alg != JetAlgorithm.AntiKt
                test_exclusive_jets_njets(cluster_seq_ptrs, cluster_seqs; njets = 4)
                test_exclusive_jets_dcut(cluster_seq_ptrs, cluster_seqs; dcut = 0.99)
            end

            for ptr in cluster_seq_ptrs
                C_JetReconstruction.jetreconstruction_ClusterSequence_free_members_(ptr)
                Libc.free(ptr)
            end
        end
    end
end
