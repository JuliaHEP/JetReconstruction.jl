"""Benchmark suite for JetReconstruction.jl compatible with PkgBenchmark.jl."""

using BenchmarkTools
using JetReconstruction

const events_file_pp = joinpath(@__DIR__, "..", "test", "data", "events.pp13TeV.hepmc3.zst")
const events_file_ee = joinpath(@__DIR__, "..", "test", "data", "events.eeH.hepmc3.zst")

const pp_events = JetReconstruction.read_final_state_particles(events_file_pp, PseudoJet)
const ee_events = JetReconstruction.read_final_state_particles(events_file_ee, EEJet)

function jet_reconstruct_harness(events; algorithm, strategy, distance, power = nothing,
                                 recombine = RecombinationMethods[RecombinationScheme.EScheme],
                                 ptmin::Real = 5.0, dcut = nothing, njets = nothing,)
    for event in events
        cs = jet_reconstruct(event; R = distance, p = power, algorithm = algorithm,
                             strategy = strategy, recombine...)
        if !isnothing(njets)
            finaljets = exclusive_jets(cs; njets = njets)
        elseif !isnothing(dcut)
            finaljets = exclusive_jets(cs; dcut = dcut)
        else
            finaljets = inclusive_jets(cs; ptmin = ptmin)
        end
    end
end

const SUITE = BenchmarkGroup()

## pp events
let pp_group = SUITE["pp"] = BenchmarkGroup(["pp"])
    for stg in [RecoStrategy.N2Plain, RecoStrategy.N2Tiled]
        strategy_group = pp_group["$stg"] = BenchmarkGroup(["$stg"])
        for alg in [JetAlgorithm.AntiKt, JetAlgorithm.CA, JetAlgorithm.Kt]
            alg_group = strategy_group["$alg"] = BenchmarkGroup(["$alg"])
            for distance in [0.4, 1.0, 1.5]
                alg_group["$distance"] = @benchmarkable jet_reconstruct_harness($pp_events;
                                                                                algorithm = $alg,
                                                                                strategy = $stg,
                                                                                distance = $distance,
                                                                                ptmin = 5.0) evals=1 samples=32
            end
        end
    end
end

## ee events
let ee_group = SUITE["ee"] = BenchmarkGroup(["ee"])
    for alg in [JetAlgorithm.Durham]
        alg_group = SUITE["ee"]["$alg"] = BenchmarkGroup(["$alg"])
        for distance in [0.4, 1.0, 1.5]
            alg_group["$distance"] = @benchmarkable jet_reconstruct_harness($ee_events;
                                                                            algorithm = $alg,
                                                                            strategy = $RecoStrategy.Best,
                                                                            distance = $distance,
                                                                            ptmin = 5.0) evals=1 samples=32
        end
    end
end
