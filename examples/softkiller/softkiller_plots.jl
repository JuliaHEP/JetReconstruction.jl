_ = """
run it with
julia --project=examples examples/softk.jl  --maxevents=100 --grid-size=0.4  --algorithm=Kt/
    --pileup-file=test/data/sk_example_PU.hepmc --hard-file=test/data/sk_example_HS.hepmc 
----------------------------------------------------------------------
"""

_ = """
A simple example of SoftKiller that reads from two HepMC3 files 
and displays plots of clustering without SoftKiller and with SoftKiller 
"""
using ArgParse
using Profile
using Logging
using JSON

using LorentzVectorHEP
using JetReconstruction
using CairoMakie

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

function process_event(event::Vector{PseudoJet}, args::Dict{Symbol, Any},
                       rapmax::Float64)
    event = select_ABS_RAP_max(event, rapmax)

    distance = args[:distance]
    algorithm = args[:algorithm]
    p = args[:power]
    strategy = args[:strategy]

    filtered_event = PseudoJet[]
    for (i, pseudo_jet) in enumerate(event)
        new_pseudo_jet = PseudoJet(JetReconstruction.px(pseudo_jet),
                                   JetReconstruction.py(pseudo_jet),
                                   JetReconstruction.pz(pseudo_jet),
                                   JetReconstruction.energy(pseudo_jet);
                                   cluster_hist_index = i)
        push!(filtered_event, new_pseudo_jet)
    end

    # Clustering of the given vector of PseudoJets 
    cluster_seq_pu = jet_reconstruct(filtered_event,
                                     R = distance, p = p, algorithm = algorithm,
                                     strategy = strategy)

    finaljets_pu_lorv = inclusive_jets(cluster_seq_pu, ptmin = 25.0)
    return finaljets_pu_lorv
end

function push_data!(event::AbstractVector, y::Vector{Float64}, phi::Vector{Float64},
                    pt::Vector{Float64}, colors::Vector{String},
                    color::String, origin::Dict{PseudoJet, String})
    for jet in event
        pj = isa(jet, PseudoJet) ? jet :
             PseudoJet(JetReconstruction.px(jet), JetReconstruction.py(jet),
                       JetReconstruction.pz(jet), JetReconstruction.energy(jet))
        push!(y, JetReconstruction.rapidity(pj))
        push!(phi, JetReconstruction.phi(pj))
        push!(pt, JetReconstruction.pt2(pj))
        if haskey(origin, jet)
            push!(colors, origin[jet] == "hard" ? "purple" : "white")
        else
            push!(colors, color)
        end
    end
end

function plot_set_up(y::Vector{Float64}, phi::Vector{Float64}, pt::Vector{Float64},
                     color::Vector{String}, plot_title::String)
    fig = Figure()
    ax = Axis(fig[1, 1],
              xlabel = "Rapidity (y)", ylabel = "Azimuthal Angle (φ)",
              title = plot_title,
              limits = ((-5.0, 5.0), (0.0, 6.5)))

    marker_sizes = 8 .+ pt ./ 50
    color_syms = Symbol.(color)

    scatter!(ax, y, phi;
             markersize = marker_sizes,
             strokecolor = :black,
             strokewidth = 1.0,
             color = color_syms,
             transparency = true)

    Legend(fig[1, 2],
           [
               PolyElement(color = :white),
               PolyElement(color = :purple),
               PolyElement(color = :royalblue3)
           ],
           ["Pileup", "Hard Event", "Jet"], "Legend")

    save_dir = "examples/softkiller"
    mkpath(save_dir)
    file_path = joinpath(save_dir, plot_title * ".png")
    save(file_path, fig)
end

function main()
    args = parse_command_line(ARGS)
    logger = ConsoleLogger(stdout, Logging.Info)
    global_logger(logger)

    if JetReconstruction.is_ee(args[:algorithm])
        jet_type = EEjet
    else
        jet_type = PseudoJet
    end
    @assert jet_type==PseudoJet "SoftKiller only supports PseudoJet"

    # Reading from input files 
    events = read_final_state_particles(args[:pileup_file],
                                        maxevents = args[:maxevents],
                                        skipevents = args[:skip],
                                        T = jet_type)

    h_events = read_final_state_particles(args[:hard_file],
                                          maxevents = args[:maxevents],
                                          skipevents = args[:skip],
                                          T = jet_type)

    # Setting dimensions for Softkiller                                   
    rapmax = 5.0
    grid_size = args[:grid_size]
    soft_killer = SoftKiller(rapmax, grid_size)

    algorithm = args[:algorithm]
    p = args[:power]

    (p, algorithm) = JetReconstruction.get_algorithm_power_consistency(p = p,
                                                                       algorithm = algorithm)
    @info "Jet reconstruction will use $(algorithm) with power $(p)"

    # This is a vector of PseudoJets that contains hard event and pileup 
    # this vector will get clustered without SoftKiller applied 
    all_jets = PseudoJet[]
    # This is a vector of PseudoJets that contains hard event and pileup 
    # but it will get clustered after SoftKiller was applied 
    all_jets_sk = PseudoJet[]

    hard_only = PseudoJet[]
    origin = Dict{PseudoJet, String}()

    # Clustering and updating data for the graphs plus filling all_jets_sk for the pile up 
    for event in events
        for pseudo_jet in event
            push!(all_jets_sk, pseudo_jet)
            push!(all_jets, pseudo_jet)
            origin[pseudo_jet] = "pileup"
        end
    end

    for pseudo_jet in h_events[2]
        push!(hard_only, pseudo_jet)
        push!(all_jets_sk, pseudo_jet)
        push!(all_jets, pseudo_jet)
        origin[pseudo_jet] = "hard"
    end

    y_hard, phi_hard, pt_hard, color_hard = Float64[], Float64[], Float64[], String[]
    push_data!(hard_only, y_hard, phi_hard, pt_hard, color_hard, "royalblue3", origin)
    plot_set_up(y_hard, phi_hard, pt_hard, color_hard,
                "Hard event only")

    y_exp, ph_exp, pt_exp, color_exp = Float64[], Float64[], Float64[], String[]
    expected = process_event(hard_only, args, rapmax)
    push_data!(expected, y_exp, ph_exp, pt_exp, color_exp, "royalblue3", origin)
    plot_set_up(y_exp, ph_exp, pt_exp, color_exp,
                "Jets, expected results")

    y_all_nosk, phi_all_nosk, pt_all_nosk, colors_all_nosk = Float64[], Float64[],
                                                             Float64[], String[]
    push_data!(all_jets, y_all_nosk, phi_all_nosk, pt_all_nosk, colors_all_nosk,
               "royalblue3", origin)
    plot_set_up(y_all_nosk, phi_all_nosk, pt_all_nosk, colors_all_nosk,
                "All PseudoJets before clustering, no SoftKiller")

    y_cl_nosk, phi_cl_nosk, pt_cl_nosk, colors_cl_nosk = Float64[], Float64[], Float64[],
                                                         String[]
    clustered_jets = process_event(all_jets, args, rapmax)
    push_data!(clustered_jets, y_cl_nosk, phi_cl_nosk, pt_cl_nosk, colors_cl_nosk,
               "royalblue3", origin)
    plot_set_up(y_cl_nosk, phi_cl_nosk, pt_cl_nosk, colors_cl_nosk,
                "All PseudoJets after clustering, no SoftKiller")

    pt_threshold = 0.00
    # Applying SoftKiller to a non-clustered vector of PseudoJets 
    reduced_event, pt_threshold = softkiller!(soft_killer, all_jets_sk)

    y_all_sk, phi_all_sk, pt_all_sk, colors_all_sk = Float64[], Float64[], Float64[],
                                                     String[]
    push_data!(reduced_event, y_all_sk, phi_all_sk, pt_all_sk, colors_all_sk, "royalblue3",
               origin)
    plot_set_up(y_all_sk, phi_all_sk, pt_all_sk, colors_all_sk,
                "All PseudoJets before clustering, with SoftKiller")

    # Clustering after Softkiller using reduced_event
    y_cl_sk, phi_cl_sk, pt_cl_sk, colors_cl_sk = Float64[], Float64[], Float64[], String[]
    clustered_jets_sk = process_event(reduced_event, args, rapmax)
    push_data!(clustered_jets_sk, y_cl_sk, phi_cl_sk, pt_cl_sk, colors_cl_sk, "royalblue3",
               origin)
    plot_set_up(y_cl_sk, phi_cl_sk, pt_cl_sk, colors_cl_sk,
                "All PseudoJets after clustering, with SoftKiller")
end

main()
