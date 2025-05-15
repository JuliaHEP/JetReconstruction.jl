#! /usr/bin/env julia
"""Use the visualisation tools to plot the reconstructed jets"""

using ArgParse
using Logging
using CairoMakie

using JetReconstruction

# Parsing for algorithm and strategy enums
include(joinpath(@__DIR__, "..", "parse-options.jl"))

function main()
    s = ArgParseSettings(autofix_names = true)
    @add_arg_table s begin
        "--event", "-e"
        help = "Event in file to visualise"
        arg_type = Int
        default = 1

        "--ptmin"
        help = "Minimum p_t for final inclusive jets (energy unit is the same as the input clusters, usually GeV)"
        arg_type = Float64
        default = 5.0

        "--exclusive-dcut"
        help = "Return all exclusive jets where further merging would have d>d_cut"
        arg_type = Float64

        "--exclusive-njets"
        help = "Return all exclusive jets once clusterisation has produced n jets"
        arg_type = Int

        "--distance", "-R"
        help = "Distance parameter for jet merging"
        arg_type = Float64
        default = 0.4

        "--algorithm", "-A"
        help = """Algorithm to use for jet reconstruction: $(join(JetReconstruction.AllJetRecoAlgorithms, ", "))"""
        arg_type = JetAlgorithm.Algorithm

        "--power", "-p"
        help = """Power value for jet reconstruction"""
        arg_type = Float64

        "--strategy", "-S"
        help = """Strategy for the algorithm, valid values: $(join(JetReconstruction.AllJetRecoStrategies, ", "))"""
        arg_type = RecoStrategy.Strategy
        default = RecoStrategy.Best

        "file"
        help = "HepMC3 event file in HepMC3 to read"
        required = true

        "output"
        help = "File for output image"
        default = "jetvis.png"
    end
    args = parse_args(ARGS, s; as_symbols = true)

    logger = ConsoleLogger(stdout, Logging.Info)
    global_logger(logger)

    events::Vector{Vector{PseudoJet}} = read_final_state_particles(args[:file],
                                                                   maxevents = args[:event],
                                                                   skipevents = args[:event])

    (p, algorithm) = JetReconstruction.get_algorithm_power_consistency(p = args[:power],
                                                                algorithm = args[:algorithm])
    cs = jet_reconstruct(events[1], R = args[:distance], p = p, algorithm = algorithm,
                         strategy = args[:strategy])

    plt = jetsplot(events[1], cs; Module = CairoMakie)
    save(args[:output], plt)

    @info "Saved jet visualisation to $(args[:output])"
end

main()
