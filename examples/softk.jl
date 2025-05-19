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
                       rapmax::Float64)
    event = select_ABS_RAP_max(event, rapmax)

    distance = args[:distance]
    algorithm = args[:algorithm]
    p = args[:power]
    strategy = args[:strategy]

    #Clustering of the given vector of PseudoJets 
    cluster_seq_pu = jet_reconstruct(event,
                                     R = distance, p = p, algorithm = algorithm,
                                     strategy = strategy)
                                     
    ize = length(event)
    println("Cluster size $(ize) \n")
    finaljets_pu_lorv = inclusive_jets(cluster_seq_pu, ptmin = 25.0)
    size = length(finaljets_pu_lorv)
    println("final size $(size) \n")
    return finaljets_pu_lorv
end

function push_data(event::AbstractVector, Ys::Vector{Float64}, Phis::Vector{Float64},
    pts::Vector{Float64}, colors::Vector{String},
    color::String, origin::Dict{PseudoJet, String})

    for (jet_i, jet) in enumerate(event)
        pj = isa(jet, PseudoJet) ? jet : PseudoJet(px(jet), py(jet), pz(jet), energy(jet))
        push!(Ys, JetReconstruction.rapidity(pj))
        push!(Phis, JetReconstruction.phi(pj))
        push!(pts, JetReconstruction.pt2(pj))
        if haskey(origin, jet)
            push!(colors, origin[jet] == "hard" ? "purple" : "white")
        else
            push!(colors, color)  # fallback if not found
        end
    end
end


function plot_set_up(Y::Vector{Float64}, Phi::Vector{Float64}, pt::Vector{Float64},
                     color::Vector{String}, plot_title::String)
    y_min, y_max = -5.0, 5.0
    phi_min, phi_max = 0.0, 6.5
    x = y_min:0.4:y_max
    y = phi_min:0.4:phi_max

    marker_sizes = 3 .+ pt/50

    format(lines) = begin
        return [isinteger(line) ? string(line) : "" for line in lines]
    end

    x_labels = format(x)
    y_labels = format(y)

    p = scatter(Y, Phi,
                xlabel = "Rapidity (y)",
                ylabel = "Azimuthal Angle (φ)",
                title = plot_title,
                markersize = marker_sizes,
                xticks = (x, x_labels),
                yticks = (y, y_labels),
                xlims = (y_min, y_max),
                ylims = (phi_min, phi_max),
                grid = true,
                framestyle = :box,
                color = color;
                alpha = 0.6,
                legend = false)

    scatter!(p, [NaN], [NaN], color = "white", markerstrokecolor = :black, label = "Pileup")
    scatter!(p, [NaN], [NaN], color = "purple", label = "Hard Event")
    scatter!(p, [NaN], [NaN], color = "royalblue3", label = "Jet")

    save_dir = "examples"
    mkpath(save_dir)
    file_path = joinpath(save_dir, plot_title * ".png")
    savefig(p, file_path)
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

    #Setting dimensions for Softkiller                                   
    rapmax = 5.0
    grid_size = args[:grid_size]
    soft_killer = SoftKiller(rapmax, grid_size)

    algorithm = args[:algorithm]
    p = args[:power]

    (p, algorithm) = JetReconstruction.get_algorithm_power_consistency(p = p,
                                                                       algorithm = algorithm)
    @info "Jet reconstruction will use $(algorithm) with power $(p)"

    #This is a vector of PseudoJets that contains hard event and pileup 
    #this vector will get clustered without SoftKiller applied 
    all_jets = PseudoJet[]
    #This is a vector of PseudoJets that contains hard event and pileup 
    #but it will get clustered after SoftKiller was applied 
    all_jets_sk = PseudoJet[]

    hard_only = PseudoJet[]
    origin= Dict{PseudoJet, String}()

    Yi, Phii, pti, colorsr = Float64[], Float64[], Float64[], String[]
    #Clustering and updating data for the graphs plus filling all_jets_sk for the pile up 
    for (ievn, event) in enumerate(events)
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
    Yh1, Phih1, pth1, colorh1 = Float64[], Float64[], Float64[], String[]
    push_data(hard_only,   Yh1, Phih1, pth1, colorh1, "royalblue3",origin)
    plot_set_up(  Yh1, Phih1, pth1, colorh1,
                "Hard event only")

    Yh, Phih, pth, colorh = Float64[], Float64[], Float64[], String[]
    expected = process_event(hard_only, args, rapmax)
    push_data(expected,  Yh, Phih, pth, colorh, "royalblue3",origin)
    plot_set_up( Yh, Phih, pth, colorh,
                "Jets, expected results")
    println("Yh: ", Yh)
    println("Phih: ", Phih)
    println("pth: ", pth)

    
    push_data(all_jets,  Yi, Phii, pti, colorsr, "royalblue3",origin)
    plot_set_up(Yi, Phii, pti, colorsr,
                "All PseudoJets before clustering, no SoftKiller")
    
    Ys, Phis, pts, colors = Float64[], Float64[], Float64[], String[]
    clusterd_jets = process_event(all_jets, args, rapmax)
    push_data(clusterd_jets,  Ys, Phis, pts, colors, "royalblue3",origin)
    plot_set_up(Ys, Phis, pts, colors,
                "All PseudoJets after clustering, no SoftKiller")


    pt_threshold = 0.00
    soft_killer_event = PseudoJet[]
    #Applying SoftKiller to a non-clustered vector of PseudoJets 
    reduced_event, pt_threshold = softkiller_apply(soft_killer, all_jets_sk, soft_killer_event, pt_threshold)

    Y_sk, Phi_sk, pt_sk, colors_sk = Float64[], Float64[], Float64[], String[]
    push_data(reduced_event,  Y_sk, Phi_sk, pt_sk, colors_sk, "royalblue3",origin)
    plot_set_up(Y_sk, Phi_sk, pt_sk, colors_sk,
                "All PseudoJets before clustering, with SoftKiller")
    
    #Clustering after Softkiller using reduced_event
    # process_event(reduced_event, args, all_jets_sk, rapmax, Ys_reduced, Phis_reduced,
    #               pts_reduced, colors_reduced, "royalblue3", origin)

    Yi3, Phii3, pti3, colors_reduced = Float64[], Float64[], Float64[], String[]
    clustered_jets_sk = process_event(soft_killer_event, args, rapmax)
    push_data(clustered_jets_sk, Yi3, Phii3, pti3, colors_reduced, "royalblue3",origin)
    println("Yi3: ", Yi3)
    println("Phii3: ", Phii3)
    println("pti3: ", pti3)
    plot_set_up(Yi3, Phii3, pti3, colors_reduced,
    "All PseudoJets after clustering, with SoftKiller")

end

main()
