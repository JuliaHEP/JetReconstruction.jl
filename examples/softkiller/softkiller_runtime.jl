using ArgParse
using Profile
using Logging
using JSON

using LorentzVectorHEP
using JetReconstruction
using Profile

# Parsing for algorithm and strategy enums
include(joinpath(@__DIR__, "..", "parse-options.jl"))

function parse_command_line(args)
    s = ArgParseSettings(autofix_names = true)
    @add_arg_table! s begin
        "--maxevents", "-n"
        help = "Maximum number of events to read. -1 to read all events from the  file."
        arg_type = Int
        default = -1

        "--skip", "-s"
        help = "Number of events to skip at beginning of the file."
        arg_type = Int
        default = 0

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

        "--dump"
        help = "Write list of reconstructed jets to a JSON formatted file"

        "--grid-size"
        help = "Size of Rectangular grid"
        arg_type = Float64
        default = 0.4

        "--hard-file"
        help = "HepMC3 event file in HepMC3 to read."
        required = true

        "--pileup-file"
        help = "HepMC3 event file in HepMC3 to read."
        required = true
    end
    return parse_args(args, s; as_symbols = true)
end

function main()
    args = parse_command_line(ARGS)
    logger = ConsoleLogger(stdout, Logging.Info)
    global_logger(logger)

    # Only PseudoJet is supported for SoftKiller
    if JetReconstruction.is_ee(args[:algorithm])
        jet_type = EEjet
    else
        jet_type = PseudoJet
    end
    @assert jet_type==PseudoJet "SoftKiller only supports PseudoJet"

    events = Vector{PseudoJet}[]

    # Reading pileup and hard event files
    events = read_final_state_particles(args[:pileup_file],
                                        maxevents = args[:maxevents],
                                        skipevents = args[:skip],
                                        T = jet_type)

    h_events = read_final_state_particles(args[:hard_file],
                                          maxevents = args[:maxevents],
                                          skipevents = args[:skip],
                                          T = jet_type)

    # Set up SoftKiller grid and rapidity range
    rapmax = 5.0
    grid_size = args[:grid_size]
    soft_killer = SoftKiller(rapmax, grid_size)

    algorithm = args[:algorithm]
    p = args[:power]

    # Ensure algorithm and power are consistent
    (p,
     algorithm) = JetReconstruction.get_algorithm_power_consistency(p = p,
                                                                    algorithm = algorithm)
    @info "Jet reconstruction will use $(algorithm) with power $(p)"

    # all_jets_sk: all PseudoJets (hard + pileup), for SoftKiller application
    all_jets_sk = PseudoJet[]

    # Fill pileup jets
    for event in events
        for pseudo_jet in event
            push!(all_jets_sk, pseudo_jet)
        end
    end

    # Fill hard event jets (from second event in h_events)
    for pseudo_jet in h_events[2]
        push!(all_jets_sk, pseudo_jet)
    end

    # Apply SoftKiller to all_jets_sk (hard + pileup) and 
    # measure the runtime of the softkiller method 
    @time softkiller!(soft_killer, all_jets_sk)
end

main()
