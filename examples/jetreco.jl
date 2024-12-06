#! /usr/bin/env julia
"""
Simple example of a jet reconstruction code that reads in a text HepMC3 file
and performs standard jet reconstruction on the final state particles.
"""

using ArgParse
using Profile
using Logging
using JSON

using LorentzVectorHEP
using JetReconstruction

# Parsing for algorithm and strategy enums
include(joinpath(@__DIR__, "parse-options.jl"))

"""
Top level call function for demonstrating the use of jet reconstruction

This uses the "jet_reconstruct" wrapper, so algorithm switching
happens inside the JetReconstruction package itself.

Final jets can be serialised if the "dump" option is given
"""
function jet_process(events::Vector{Vector{T}};
                     distance::Real = 0.4,
                     algorithm::Union{JetAlgorithm.Algorithm, Nothing} = nothing,
                     p::Union{Real, Nothing} = nothing,
                     ptmin::Real = 5.0,
                     dcut = nothing,
                     njets = nothing,
                     strategy::RecoStrategy.Strategy,
                     dump::Union{String, Nothing} = nothing) where {T <:
                                                                    JetReconstruction.FourMomentum}

    # Set consistent algorithm and power
    (p, algorithm) = JetReconstruction.get_algorithm_power_consistency(p = p,
                                                                       algorithm = algorithm)
    @info "Jet reconstruction will use $(algorithm) with power $(p)"

    # A friendly label for the algorithm and final jet selection
    if !isnothing(njets)
        @info "Will produce exclusive jets with n_jets = $(njets)"
    elseif !isnothing(dcut)
        @info "Will produce exclusive jets with dcut = $(dcut)"
    else
        @info "Will produce inclusive jets with ptmin = $(ptmin)"
    end

    # Now run over each event
    for (ievt, event) in enumerate(events)
        # Run the jet reconstruction
        cluster_seq = jet_reconstruct(event, R = distance, p = p, algorithm = algorithm,
                                      strategy = strategy)
        # Now select jets, with inclusive or exclusive parameters
        if !isnothing(njets)
            finaljets = exclusive_jets(cluster_seq; njets = njets)
        elseif !isnothing(dcut)
            finaljets = exclusive_jets(cluster_seq; dcut = dcut)
        else
            finaljets = inclusive_jets(cluster_seq; ptmin = ptmin)
        end
        @info begin
            jet_output = "Event $(ievt)\n"
            sort!(finaljets, by = x -> pt(x), rev = true)
            for (ijet, jet) in enumerate(finaljets)
                jet_output *= " $(ijet) - $(jet)\n"
            end
            "$(jet_output)"
        end
        if !isnothing(dump)
            push!(jet_collection, FinalJets(ievt, finaljets))
        end
    end

    if !isnothing(dump)
        open(dump, "w") do io
            JSON.print(io, jet_collection, 2)
        end
    end
end

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

        "file"
        help = "HepMC3 event file in HepMC3 to read."
        required = true
    end
    return parse_args(args, s; as_symbols = true)
end

function main()
    args = parse_command_line(ARGS)
    logger = ConsoleLogger(stdout, Logging.Info)
    global_logger(logger)
    # Try to read events into the correct type!
    # If we don't have an algorithm we default to PseudoJet
    if !isnothing(args[:algorithm])
        JetReconstruction.is_ee(args[:algorithm])
        jet_type = EEjet{Float64}
    else
        jet_type = PseudoJet{Float64}
    end
    events::Vector{Vector{jet_type}} = read_final_state_particles(args[:file],
                                                                  maxevents = args[:maxevents],
                                                                  skipevents = args[:skip],
                                                                  T = jet_type)
    if isnothing(args[:algorithm]) && isnothing(args[:power])
        @warn "Neither algorithm nor power specified, defaulting to pp event AntiKt"
        args[:algorithm] = JetAlgorithm.AntiKt
    end
    jet_process(events, distance = args[:distance], algorithm = args[:algorithm],
                p = args[:power],
                strategy = args[:strategy],
                ptmin = args[:ptmin], dcut = args[:exclusive_dcut],
                njets = args[:exclusive_njets],
                dump = args[:dump])
    nothing
end

main()
