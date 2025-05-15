#! /usr/bin/env julia
"""Use the visualisation tools to produce an animation of the jet reconstruction"""

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
        help = "File for output animation"
        default = "jetreco.mp4"
    end
    args = parse_args(ARGS, s; as_symbols = true)

    logger = ConsoleLogger(stdout, Logging.Info)
    global_logger(logger)

    events::Vector{Vector{PseudoJet}} = read_final_state_particles(args[:file],
                                                                   maxevents = args[:event],
                                                                   skipevents = args[:event])
    if isnothing(args[:algorithm]) && isnothing(args[:power])
        @warn "Neither algorithm nor power specified, defaulting to AntiKt"
        args[:algorithm] = JetAlgorithm.AntiKt
    end
    # Set consistent algorithm and power
    (p, algorithm) = JetReconstruction.get_algorithm_power_consistency(p = args[:power],
                                                                algorithm = args[:algorithm])
    cs = jet_reconstruct(events[1], R = args[:distance], p = p, algorithm = algorithm,
                         strategy = args[:strategy])

    animatereco(cs, args[:output]; azimuth = (1.8, 3.0), elevation = 0.5,
                perspective = 0.5, framerate = 20, ancestors = true,
                barsize_phi = 0.1, barsize_y = 0.3)

    @info "Saved jet reconstruction animation to $(args[:output])"
end

main()
