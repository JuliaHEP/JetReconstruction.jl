#! /usr/bin/env julia
"""
Wrapper to run jet reco code feeding in the standard set of HepMC events that
are used for benchmarking jet reco.
"""

using FlameGraphs: FlameGraphs
using ProfileSVG: ProfileSVG

using ArgParse
using Profile
using Colors
using StatProfilerHTML
using Logging
using JSON

using LorentzVectorHEP
using JetReconstruction

function profile_code(jet_reconstruction, events, niters)
	Profile.init(n = 5*10^6, delay = 0.00001)
	profile_events(events) = begin
		for evt in events
			jet_reconstruction(evt, R = 0.4)
		end
	end
	profile_events(events[1:1])
	@profile for i ∈ 1:niters
		profile_events(events)
	end
	statprofilehtml()
	fcolor = FlameGraphs.FlameColors(
		reverse(colormap("Blues", 15))[1:5],
		colorant"slategray4",
		colorant"gray95",
		reverse(colormap("Reds", 15))[1:5],
		reverse(sequential_palette(39, 10; s = 38, b = 2))[1:5],#brownish pallette
	)
	ProfileSVG.save(
		fcolor,
		joinpath("statprof", "profsvg.svg");
		combine = true,
		timeunit = :ms,
		font = "Arial, Helvetica, sans-serif",
	)
	println(
		"Flame graph from ProfileSVG.jl at file://",
		abspath("statprof/profsvg.svg"),
		"\n",
		"""
\tRed tint:          Runtime dispatch
\tBrown/yellow tint: Garbage collection
\tBlue tint:         OK
""",
	)
end

function jet_process(
	events::Vector{Vector{PseudoJet}};
	ptmin::Float64 = 5.0,
	distance::Float64 = 0.4,
	power::Integer = -1,
	strategy::JetRecoStrategy,
	nsamples::Integer = 1,
	gcoff::Bool = false,
	profile::Bool = false,
    alloc::Bool = false,
	dump::Union{String, Nothing} = nothing,
)
	@info "Will process $(size(events)[1]) events"

	# Strategy
	if (strategy == N2Plain)
		jet_reconstruction = plain_jet_reconstruct
	elseif (strategy == N2Tiled || stragegy == Best)
		jet_reconstruction = tiled_jet_reconstruct
	else
		throw(ErrorException("Strategy not yet implemented"))
	end

	# Build internal EDM structures for timing measurements, if needed
	# For N2Tiled and N2Plain this is unnecessary as both these reconstruction
	# methods can process PseudoJets directly
	if (strategy == N2Tiled) || (strategy == N2Plain)
		event_vector = events
	end

	# If we are dumping the results, setup the JSON structure
	if !isnothing(dump)
		jet_collection = FinalJets[]
	end

	# Warmup code if we are doing a multi-sample timing run
	if nsamples > 1 || profile
		@info "Doing initial warm-up run"
		for event in event_vector
			finaljets, _ = jet_reconstruction(event, R = distance, p = power)
			final_jets(finaljets, ptmin)
		end
	end

	if profile
		profile_code(jet_reconstruction, event_vector, nsamples)
		return nothing
	end

    if alloc
        println("Memory allocation statistics:")
        @timev for event in event_vector
            finaljets, _ = jet_reconstruction(event, R = distance, p = power)
			final_jets(finaljets, ptmin)
        end
        return nothing
    end

	# Now setup timers and run the loop
	GC.gc()
	cummulative_time = 0.0
	cummulative_time2 = 0.0
	lowest_time = typemax(Float64)
	for irun ∈ 1:nsamples
		gcoff && GC.enable(false)
		t_start = time_ns()
		for (ievt, event) in enumerate(event_vector)
			finaljets, _ = jet_reconstruction(event, R = distance, p = power, ptmin=ptmin)
			fj = final_jets(finaljets, ptmin)
			# Only print the jet content once
			if irun == 1
				@info begin
					jet_output = "Event $(ievt)\n"
					for (ijet, jet) in enumerate(fj)
						jet_output *= " $(ijet) - $(jet)\n"
					end
					"$(jet_output)"
				end
				if !isnothing(dump)
					push!(jet_collection, FinalJets(ievt, fj))
				end
			end
		end
		t_stop = time_ns()
		gcoff && GC.enable(true)
		dt_μs = convert(Float64, t_stop - t_start) * 1.e-3
        if nsamples > 1
			@info "$(irun)/$(nsamples) $(dt_μs)"
		end
		cummulative_time += dt_μs
		cummulative_time2 += dt_μs^2
		lowest_time = dt_μs < lowest_time ? dt_μs : lowest_time
	end

	mean = cummulative_time / nsamples
	cummulative_time2 /= nsamples
	if nsamples > 1
		sigma = sqrt(nsamples / (nsamples - 1) * (cummulative_time2 - mean^2))
	else
		sigma = 0.0
	end
	mean /= length(events)
	sigma /= length(events)
	lowest_time /= length(events)
	# Why also record the lowest time? 
	# 
	# The argument is that on a "busy" machine, the run time of an application is
	# always TrueRunTime+Overheads, where Overheads is a nuisance parameter that
	# adds jitter, depending on the other things the machine is doing. Therefore
	# the minimum value is (a) more stable and (b) reflects better the intrinsic
	# code performance.
	println("Processed $(length(events)) events $(nsamples) times")
	println("Average time per event $(mean) ± $(sigma) μs")
	println("Lowest time per event $lowest_time μs")

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
		help = "Minimum p_t for final jets (GeV)"
		arg_type = Float64
		default = 5.0

		"--distance", "-R"
		help = "Distance parameter for jet merging"
		arg_type = Float64
		default = 0.4

		"--power"
		help = "Distance measure momentum power (-1 - antikt; 0 - Cambridge/Achen; 1 - inclusive k_t)"
		arg_type = Int
		default = -1

		"--strategy"
		help = "Strategy for the algorithm, valid values: Best, N2Plain, N2Tiled, N2TiledSoAGlobal, N2TiledSoATile"
		arg_type = JetRecoStrategy
		default = N2Plain

		"--nsamples", "-m"
		help = "Number of measurement points to acquire."
		arg_type = Int
		default = 1

		"--gcoff"
		help = "Turn off Julia garbage collector during each time measurement."
		action = :store_true

		"--profile"
		help = "Profile code and generate a flame graph."
		action = :store_true

        "--alloc"
        help = "Provide memory allocation statistics."
        action = :store_true

		"--dump"
		help = "Write list of recontructed jets to a JSON formatted file"

		"--info"
		help = "Print info level log messages"
		action = :store_true

		"--debug"
		help = "Print debug level log messages"
		action = :store_true

		"file"
		help = "HepMC3 event file in HepMC3 to read."
		required = true
	end
	return parse_args(args, s; as_symbols = true)
end


function ArgParse.parse_item(::Type{JetRecoStrategy}, x::AbstractString)
	if (x == "Best")
		return JetRecoStrategy(0)
	elseif (x == "N2Plain")
		return JetRecoStrategy(1)
	elseif (x == "N2Tiled")
		return JetRecoStrategy(2)
	else
		throw(ErrorException("Invalid value for strategy: $(x)"))
	end
end

main() = begin
	args = parse_command_line(ARGS)
	if args[:debug]
		logger = ConsoleLogger(stdout, Logging.Debug)
	elseif args[:info]
		logger = ConsoleLogger(stdout, Logging.Info)
	else
		logger = ConsoleLogger(stdout, Logging.Warn)
	end
	global_logger(logger)
	events::Vector{Vector{PseudoJet}} =
		read_final_state_particles(args[:file], maxevents = args[:maxevents], skipevents = args[:skip])
	jet_process(events, ptmin = args[:ptmin], distance = args[:distance], 
        power = args[:power], strategy = args[:strategy],
		nsamples = args[:nsamples], gcoff = args[:gcoff], profile = args[:profile],
		alloc = args[:alloc], dump = args[:dump])
	nothing
end

# Running through the debugger in VS Code actually has
# ARGS[0] = "/some/path/.vscode/extensions/julialang.language-julia-1.47.2/scripts/debugger/run_debugger.jl",
# so then the program does nothing at all if it only tests for abspath(PROGRAM_FILE) == @__FILE__
# In addition, deep debugging with Infiltrator needs to start in an interactive session,
# so execute main() unconditionally
main()
