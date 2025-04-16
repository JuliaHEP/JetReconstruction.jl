
using ArgParse
using Profile
using Logging
using JSON

using LorentzVectorHEP
using JetReconstruction
using Plots 
using Pkg
Pkg.develop(path="/Users/emadimtrova/Desktop/uni/spring25/UTRA/JetReconstruction.jl")

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


function process_event(event::Vector{PseudoJet}, args::Dict{Symbol, Any}, all_jets::Vector{PseudoJet},
    rapmax::Float64, Ys::Vector{Float64}, Phis::Vector{Float64},
    pts::Vector{Float64}, colors::Vector{String}, 
    color::String, jet_origin::Vector{String})

    event = select_ABS_RAP_max(event, rapmax)

    distance = args[:distance] 
    algorithm = args[:algorithm]
    p = args[:power]
    strategy = args[:strategy]

    cluster_seq_pu = jet_reconstruct(event, 
                                     R = distance, p = p, algorithm = algorithm,
                                    strategy = strategy)
    
    #clustering 
    finaljets_pu_lorv = inclusive_jets(cluster_seq_pu)
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


main() = begin 
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

    events = read_final_state_particles(args[:pileup_file],
                                        maxevents = args[:maxevents],
                                        skipevents = args[:skip],
                                        T = jet_type)


    h_events = read_final_state_particles(args[:hard_file],
                                        maxevents = args[:maxevents],
                                        skipevents = args[:skip],
                                        T = jet_type)

                                      
    rapmax = 5.0 
    grid_size = args[:grid_size]
    soft_killer = SoftKiller(rapmax, grid_size)

    algorithm = args[:algorithm]
    p = args[:power]

    (p, algorithm) = JetReconstruction.get_algorithm_power_consistency(p = p,
    algorithm = algorithm)
    @info "Jet reconstruction will use $(algorithm) with power $(p)"

    all_events = PseudoJet[]
    jet_origin = String[]

    hard_event = h_events[1]
    hard_event = select_ABS_RAP_max(hard_event,rapmax)

    Ys, Phis, pts = Float64[], Float64[], Float64[]
    colors = String[]

    #all events post clutering 
    process_event(h_events[1], args, all_events, rapmax, Ys, Phis, pts, colors, "purple", jet_origin)

    n_events = length(events)
    println("EVENT: $n_events")
    for (ievn, event) in enumerate(events)
        process_event(event, args, all_events, rapmax, Ys, Phis, pts, colors, "black", jet_origin)        
    end 
    
    num_events = length(all_events)
    println("EVENT BEFORE: $num_events")

    pt_threshold = 0.00
    soft_killer_event = PseudoJet[]

    plot_set_up(Ys, Phis, pts, colors, "All Jets Before SoftKiller")

    reduced_event, pt_threshold = apply(soft_killer, all_events, soft_killer_event, pt_threshold)
    
    Ys_reduced, Phis_reduced, pts_reduced = Float64[], Float64[], Float64[]
    colors_reduced = String[]

    for jet in reduced_event
        push!(Ys_reduced, JetReconstruction.rapidity(jet))
        push!(Phis_reduced, JetReconstruction.phi(jet))
        push!(pts_reduced, JetReconstruction.pt2(jet))
    
        idx = findfirst(x -> x === jet, all_events)
        if idx === nothing
            push!(colors_reduced, "pink")
        else
            push!(colors_reduced, jet_origin[idx] == "hard" ? "purple" : "black")
        end
    end

    num_events_after = length(reduced_event)
    println("EVENT AFTER: $num_events_after")
    #process_event(reduced_event, args, all_events, rapmax, Ys_reduced, Phis_reduced, pts_reduced, colors_reduced, "green", jet_origin)
    
    plot_set_up(Ys_reduced, Phis_reduced, pts_reduced, colors_reduced, "All Jets After SoftKiller and Clustering")

end


main()