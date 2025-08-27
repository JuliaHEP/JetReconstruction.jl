using ArgParse
using Logging
using JSON

using LorentzVectorHEP
using JetReconstruction

# Parsing for algorithm and strategy enums
include(joinpath(@__DIR__, "..", "parse-options.jl"))

function parse_command_line(args)
    s = ArgParseSettings(autofix_names = true)
    @add_arg_table! s begin
        "--pileup-maxevents", "-n"
        help = "Maximum number of pileup events to read. -1 to read all events from the pileup events file."
        arg_type = Int
        default = -1

        "--pileup-skip", "-s"
        help = "Number of events to skip in the pileup events file."
        arg_type = Int
        default = 0

        "--eventno", "-e"
        help = "Event number to process from the hard scatter events file. If not specified, the first event will be processed."
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

        "--grid-size"
        help = "Size of Rectangular grid"
        arg_type = Float64
        default = 0.4

        "hard-file"
        help = "HepMC3 event file containing hard scatter events."
        required = true

        "pileup-file"
        help = "HepMC3 event file containing pileup events."
        required = true
    end
    return parse_args(args, s; as_symbols = true)
end

function main()
    args = parse_command_line(ARGS)
    logger = ConsoleLogger(stdout, Logging.Info)
    global_logger(logger)

    # Only PseudoJet is supported for SoftKiller
    @assert JetReconstruction.is_pp(args[:algorithm]) "SoftKiller only supports pp algorithms and PseudoJet"
    jet_type = PseudoJet

    events = Vector{PseudoJet}[]

    args[:pileup_file] = normpath(joinpath(@__DIR__, args[:pileup_file]))
    args[:hard_file] = normpath(joinpath(@__DIR__, args[:hard_file]))

    # Reading pileup and hard event files
    events = read_final_state_particles(args[:pileup_file], jet_type;
                                        maxevents = args[:pileup_maxevents],
                                        skipevents = args[:pileup_skip])

    h_events = read_final_state_particles(args[:hard_file], jet_type;
                                          maxevents = args[:eventno],
                                          skipevents = args[:eventno])

    # Set up SoftKiller grid and rapidity range
    rapmax = 5.0
    grid_size = args[:grid_size]
    soft_killer = SoftKiller(rapmax, grid_size)

    # Ensure algorithm and power are consistent
    if isnothing(args[:algorithm]) && isnothing(args[:power])
        @warn "Neither algorithm nor power specified, defaulting to AntiKt"
        args[:algorithm] = JetAlgorithm.AntiKt
    end

    # all_jets_sk: all PseudoJets (hard + pileup), for SoftKiller application
    all_jets_sk = PseudoJet[]

    # Fill pileup jets
    for event in events
        append!(all_jets_sk, event)
    end

    # Fill hard event jets from first read event in h_events (actually should be only one event!)
    append!(all_jets_sk, h_events[1])

    # Apply SoftKiller to all_jets_sk (hard + pileup)
    reduced_event, pt_threshold = softkiller(soft_killer, all_jets_sk)
    @info "SoftKiller applied: $(length(reduced_event)) clusters remaining from $(length(all_jets_sk)), pt threshold = $pt_threshold"

    # Workaround for cluster_hist_indedx not being correct after SoftKiller filtering
    fixed_reduced_event = PseudoJet[]
    sizehint!(fixed_reduced_event, length(reduced_event))
    for i in eachindex(reduced_event)
        push!(fixed_reduced_event,
              PseudoJet(reduced_event[i].px, reduced_event[i].py, reduced_event[i].pz,
                        reduced_event[i].E, i,
                        reduced_event[i]._pt2, reduced_event[i]._inv_pt2,
                        reduced_event[i]._rap, reduced_event[i]._phi))
    end
    cs = jet_reconstruct(fixed_reduced_event; algorithm = args[:algorithm],
                         R = args[:distance], p = args[:power],
                         strategy = args[:strategy])
    @info "Reconstructed softkiller filtered clusters with algorithm $(args[:algorithm]), radius $(args[:distance]) and strategy $(args[:strategy])"
end

main()
