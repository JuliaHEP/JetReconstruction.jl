"""
run it with
julia --project=examples examples/softk.jl  --maxevents=100 --grid-size=0.4  --algorithm=Kt/
    --pileup-file=test/data/sk_example_PU.hepmc --hard-file=test/data/sk_example_HS.hepmc 
----------------------------------------------------------------------
"""

"""
A simple example of SoftKiller that read from two HepMC3 files 
and displays plots of clustering without SoftKiller and with SoftKiller 
"""
using ArgParse
using Profile
using Logging
using JSON

using LorentzVectorHEP
using JetReconstruction
using Plots

include(joinpath(@__DIR__, "parse-options.jl"))

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

function process_event(event::Vector{PseudoJet}, args::Dict{Symbol, Any},
                       all_jets::Vector{PseudoJet},
                       rapmax::Float64, Ys::Vector{Float64}, Phis::Vector{Float64},
                       pts::Vector{Float64}, colors::Vector{String},
                       color::String, jet_origin::Vector{String})
    event = select_ABS_RAP_max(event, rapmax)

    distance = args[:distance]
    algorithm = args[:algorithm]
    p = args[:power]
    strategy = args[:strategy]

    #Clustering of the given vector of PseudoJets 
    cluster_seq_pu = jet_reconstruct(event,
                                     R = distance, p = p, algorithm = algorithm,
                                     strategy = strategy)

    finaljets_pu_lorv = inclusive_jets(cluster_seq_pu)

    #Updating information for plots
    for jet in finaljets_pu_lorv
        pj = PseudoJet(px(jet), py(jet), pz(jet), energy(jet))
        push!(all_jets, pj)
        push!(jet_origin, color == "purple" ? "hard" : "pileup")
        push!(Ys, JetReconstruction.rapidity(pj))
        push!(Phis, JetReconstruction.phi(pj))
        push!(pts, JetReconstruction.pt2(pj))
        push!(colors, color)
    end
end

function main()
    `read in input particles
    //
    // since we use here simulated data we can split the hard event
    // from the full (i.e. with pileup added) one
    //`
    args = parse_command_line(ARGS)
    logger = ConsoleLogger(stdout, Logging.Info)
    global_logger(logger)

    if JetReconstruction.is_ee(args[:algorithm])
        jet_type = EEjet
    else
        jet_type = PseudoJet
    end

    hard_event = PseudoJet[]
    events = Vector{PseudoJet}[]

    #Reading from input files 
    events = read_final_state_particles(args[:pileup_file],
                                        maxevents = args[:maxevents],
                                        skipevents = args[:skip],
                                        T = jet_type)

    h_events = read_final_state_particles(args[:hard_file],
                                          maxevents = args[:maxevents],
                                          skipevents = args[:skip],
                                          T = jet_type)

    #Setting dimetions for Softkiller                                   
    rapmax = 5.0
    grid_size = args[:grid_size]
    soft_killer = SoftKiller(rapmax, grid_size)

    algorithm = args[:algorithm]
    p = args[:power]

    (p,
     algorithm) = JetReconstruction.get_algorithm_power_consistency(p = p,
                                                                    algorithm = algorithm)
    @info "Jet reconstruction will use $(algorithm) with power $(p)"

    #This is a vector of PseudoJets that contatins hard event and pileup 
    #this vector will get clustered without SoftKiller applied 
    all_jets = PseudoJet[]
    #This is a vector of PseudoJets that contatins hard event and pileup 
    #but it will get clustered after SoftKiller was applied 
    all_jets_sk = PseudoJet[]
    #This is a vecotor of type string to track if a particle came from hard_events 
    #or from events - this is to mainain accurate coloring for the "Before" graph 
    jet_origin = String[]

    hard_event = h_events[1]

    #keep the particles up to 5 units in rapidity
    hard_event = select_ABS_RAP_max(hard_event, rapmax)

    Ys, Phis, pts = Float64[], Float64[], Float64[]
    colors = String[]

    #Clustering and updating data for the gaphs 
    process_event(h_events[1], args, all_jets, rapmax, Ys, Phis, pts, colors, "purple",
                  jet_origin)

    #Filling all_jets_sk for hard event
    for pseudo_jet in h_events[1]
        push!(all_jets_sk, pseudo_jet)
    end

    n_events = length(events)
    println("EVENT: $n_events")

    #Clustering and updating data for the gaphs plus filling all_jets_sk for the pile up 
    for (ievn, event) in enumerate(events)
        for pseudo_jet in event
            push!(all_jets_sk, pseudo_jet)
        end
        process_event(event, args, all_jets, rapmax, Ys, Phis, pts, colors, "black",
                      jet_origin)
    end

    num_events = length(all_jets)
    println("EVENT BEFORE: $num_events")
    #Plotting of clustering without SoftKiller applied 
    plot_set_up(Ys, Phis, pts, colors, "All Jets Before SoftKiller")

    pt_threshold = 0.00
    soft_killer_event = PseudoJet[]

    #Applying SoftKiller to a non-clistered vector of PseudoJets 
    reduced_event,
    pt_threshold = apply(soft_killer, all_jets_sk, soft_killer_event, pt_threshold)

    Ys_reduced, Phis_reduced, pts_reduced = Float64[], Float64[], Float64[]
    colors_reduced = String[]

    num_events_after = length(reduced_event)
    println("EVENT AFTER: $num_events_after")

    #Clustering after Softkiller using reduced_event
    process_event(reduced_event, args, all_jets_sk, rapmax, Ys_reduced, Phis_reduced,
                  pts_reduced, colors_reduced, "royalblue3", jet_origin)

    #Plotting reduced_event
    plot_set_up(Ys_reduced, Phis_reduced, pts_reduced, colors_reduced,
                "All Jets After SoftKiller and Clustering")
end

main()
