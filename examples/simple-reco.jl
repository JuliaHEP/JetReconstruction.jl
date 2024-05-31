#! /usr/bin/env julia
"""
Simple example of a jet reconstruction code that reads in a text HepMC3 file
and performs standard jet reconstruction on the final state particles.
"""

using ArgParse
using Profile
using Colors
using StatProfilerHTML
using Logging
using JSON

using LorentzVectorHEP
using JetReconstruction

"""
Top level call funtion for demonstrating the use of jet reconstruction

This uses the "jet_reconstruct" wrapper, so algorithm switching
happens inside the JetReconstruction package itself.

Final jets can be serialised if the "dump" option is given
"""
function jet_process(
	events::Vector{Vector{PseudoJet}};
	distance::Real = 0.4,
	power::Integer = -1,
	ptmin::Real = 5.0,
	dcut = nothing,
	njets = nothing,
	strategy::RecoStrategy.Strategy,
	dump::Union{String, Nothing} = nothing,
)
	@info "Will process $(size(events)[1]) events"

	# If we are dumping the results, setup the JSON structure
	if !isnothing(dump)
		jet_collection = FinalJets[]
	end

    # A friendly label for the algorithm and strategy
    if !isnothing(njets)
        @info "Running exclusive jets with n_jets = $(njets)"
    elseif !isnothing(dcut)
        @info "Running exclusive jets with dcut = $(dcut)"
    else
        @info "Running inclusive jets with ptmin = $(ptmin)"
    end

	# Now run over each event
    for (ievt, event) in enumerate(events)
        if !isnothing(njets)
            finaljets = exclusive_jets(jet_reconstruct(event, R = distance, p = power, strategy = strategy), njets=njets)
        elseif !isnothing(dcut)
            finaljets = exclusive_jets(jet_reconstruct(event, R = distance, p = power, strategy = strategy), dcut=dcut)
        else
            finaljets = inclusive_jets(jet_reconstruct(event, R = distance, p = power, strategy = strategy), ptmin)
        end
        @info begin
            jet_output = "Event $(ievt)\n"
            sort!(finaljets, by = x -> pt(x), rev=true)
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

parse_command_line(args) = begin
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

		"--power"
		help = "Distance measure momentum power (-1 - antikt; 0 - Cambridge/Aachen; 1 - inclusive k_t)"
		arg_type = Int
		default = -1

		"--strategy"
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


function ArgParse.parse_item(::Type{RecoStrategy.Strategy}, x::AbstractString)
	s = tryparse(RecoStrategy.Strategy, x)
	if s === nothing
		throw(ErrorException("Invalid value for strategy: $(x)"))
	end
	s
end

main() = begin
	args = parse_command_line(ARGS)
	logger = ConsoleLogger(stdout, Logging.Info)
	global_logger(logger)
	events::Vector{Vector{PseudoJet}} =
		read_final_state_particles(args[:file], maxevents = args[:maxevents], skipevents = args[:skip])
	jet_process(events, distance = args[:distance], power = args[:power], strategy = args[:strategy],
		ptmin = args[:ptmin], dcut = args[:exclusive_dcut], njets = args[:exclusive_njets],
		dump = args[:dump])
	nothing
end

main()
