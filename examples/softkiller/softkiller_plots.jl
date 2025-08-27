using ArgParse
using Logging
using JSON

using LorentzVectorHEP
using JetReconstruction
using CairoMakie

# Parsing for algorithm and strategy enums
include(joinpath(@__DIR__, "..", "parse-options.jl"))

function parse_command_line(args)
    s = ArgParseSettings(autofix_names = true)
    @add_arg_table! s begin
        "--pileup-maxevents", "-n"
        help = "Maximum number of events to read. -1 to read all events from the pileup events file."
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
        help = "HepMC3 event file in HepMC3 to read."
        required = true

        "pileup-file"
        help = "HepMC3 event file in HepMC3 to read."
        required = true
    end
    return parse_args(args, s; as_symbols = true)
end

# Cluster a single event and return inclusive jets above ptmin
function process_event(event::Vector{PseudoJet}, args::Dict{Symbol, Any},
                       rapmax::Float64)
    event = select_ABS_RAP_max(event, rapmax)

    distance = args[:distance]
    algorithm = args[:algorithm]
    p = args[:power]
    strategy = args[:strategy]

    filtered_event = PseudoJet[]
    for (i, pseudo_jet) in enumerate(event)
        # Reconstruct PseudoJet with cluster_hist_index for tracking
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

# Helper to extract rapidity, phi, pt2, and color for plotting
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
        # Color code by origin (hard/pileup) if available
        if haskey(origin, jet)
            push!(colors, origin[jet] == "hard" ? "purple" : "white")
        else
            push!(colors, color)
        end
    end
end

# Plot a scatter of jets or particles in (y, phi) with pt-dependent marker size
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

    save(plot_title * ".png", fig)
    @info "Plot '$(plot_title)' saved"
end

function main()
    args = parse_command_line(ARGS)
    logger = ConsoleLogger(stdout, Logging.Info)
    global_logger(logger)

    # Only PseudoJet is supported for SoftKiller
    @assert JetReconstruction.is_pp(args[:algorithm]) "SoftKiller only supports pp algorithms and PseudoJet"
    jet_type = PseudoJet

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

    # all_jets: all PseudoJets (hard + pileup), before SoftKiller
    # all_jets_sk: all PseudoJets (hard + pileup), for SoftKiller application
    # hard_only: only hard event PseudoJets
    # origin: maps PseudoJet to "hard" or "pileup" for coloring
    all_jets = PseudoJet[]
    all_jets_sk = PseudoJet[]
    hard_only = PseudoJet[]
    origin = Dict{PseudoJet, String}()

    # Fill pileup jets
    for event in events
        for pseudo_jet in event
            push!(all_jets_sk, pseudo_jet)
            push!(all_jets, pseudo_jet)
            origin[pseudo_jet] = "pileup"
        end
    end

    # Fill hard event jets (only the first event will be read)
    for pseudo_jet in h_events[1]
        push!(hard_only, pseudo_jet)
        push!(all_jets_sk, pseudo_jet)
        push!(all_jets, pseudo_jet)
        origin[pseudo_jet] = "hard"
    end

    # Plot hard event only
    y_hard, phi_hard, pt_hard, color_hard = Float64[], Float64[], Float64[], String[]
    push_data!(hard_only, y_hard, phi_hard, pt_hard, color_hard, "royalblue3", origin)
    plot_set_up(y_hard, phi_hard, pt_hard, color_hard,
                "Hard event only")

    # Plot expected jets from hard event
    y_exp, ph_exp, pt_exp, color_exp = Float64[], Float64[], Float64[], String[]
    expected = process_event(hard_only, args, rapmax)
    push_data!(expected, y_exp, ph_exp, pt_exp, color_exp, "royalblue3", origin)
    plot_set_up(y_exp, ph_exp, pt_exp, color_exp,
                "Jets, expected results")

    # Plot all PseudoJets before clustering, no SoftKiller
    y_all_nosk, phi_all_nosk, pt_all_nosk,
    colors_all_nosk = Float64[], Float64[],
                      Float64[], String[]
    push_data!(all_jets, y_all_nosk, phi_all_nosk, pt_all_nosk, colors_all_nosk,
               "royalblue3", origin)
    plot_set_up(y_all_nosk, phi_all_nosk, pt_all_nosk, colors_all_nosk,
                "All PseudoJets before clustering, no SoftKiller")

    # Plot all PseudoJets after clustering, no SoftKiller
    y_cl_nosk, phi_cl_nosk, pt_cl_nosk,
    colors_cl_nosk = Float64[], Float64[], Float64[],
                     String[]
    clustered_jets = process_event(all_jets, args, rapmax)
    push_data!(clustered_jets, y_cl_nosk, phi_cl_nosk, pt_cl_nosk, colors_cl_nosk,
               "royalblue3", origin)
    plot_set_up(y_cl_nosk, phi_cl_nosk, pt_cl_nosk, colors_cl_nosk,
                "All PseudoJets after clustering, no SoftKiller")

    pt_threshold = 0.00
    # Apply SoftKiller to all_jets_sk (hard + pileup)
    reduced_event, pt_threshold = softkiller(soft_killer, all_jets_sk)
    @info "SoftKiller applied: $(length(reduced_event)) clusters remaining from $(length(all_jets_sk)), pt threshold = $pt_threshold"

    # Plot all PseudoJets after SoftKiller, before clustering
    y_all_sk, phi_all_sk, pt_all_sk,
    colors_all_sk = Float64[], Float64[], Float64[],
                    String[]
    push_data!(reduced_event, y_all_sk, phi_all_sk, pt_all_sk, colors_all_sk, "royalblue3",
               origin)
    plot_set_up(y_all_sk, phi_all_sk, pt_all_sk, colors_all_sk,
                "All PseudoJets before clustering, with SoftKiller")

    # Plot all PseudoJets after SoftKiller and clustering
    y_cl_sk, phi_cl_sk, pt_cl_sk, colors_cl_sk = Float64[], Float64[], Float64[], String[]
    clustered_jets_sk = process_event(reduced_event, args, rapmax)
    push_data!(clustered_jets_sk, y_cl_sk, phi_cl_sk, pt_cl_sk, colors_cl_sk, "royalblue3",
               origin)
    plot_set_up(y_cl_sk, phi_cl_sk, pt_cl_sk, colors_cl_sk,
                "All PseudoJets after clustering, with SoftKiller")
end

main()
