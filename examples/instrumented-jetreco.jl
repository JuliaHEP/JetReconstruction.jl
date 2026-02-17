#! /usr/bin/env julia
"""
Example of running jet reconstruction on a HepMC3 file, which also
produces timing measurements.

Options are available to profile the code (flamegraphs and allocations).
"""

using FlameGraphs: FlameGraphs

using ArgParse
using Profile
using Printf
using StatProfilerHTML
using Statistics
using StatsBase
using Logging
using JSON
using UnicodePlots

using LorentzVectorHEP
using JetReconstruction

# Parsing for algorithm and strategy enums
include(joinpath(@__DIR__, "parse-options.jl"))

"""
    profile_code(events::Vector{Vector{T}}, profile_subdir, nsamples; R = 0.4, p = -1,
                      algorithm::JetAlgorithm.Algorithm = JetAlgorithm.AntiKt,
                      strategy = RecoStrategy.N2Tiled,
                      recombine = RecombinationMethods[RecombinationScheme.ESchemeRaw]) where {T <:
                                                              JetReconstruction.FourMomentum}

Profile the jet reconstruction code using the `@profile` macro and generate a
flamegraph which is saved to the `profile/profile_subdir` directory.
"""
function profile_code(events::Vector{Vector{T}}, profile, nsamples; R = 0.4, p = -1,
                      γ::Union{Real, Nothing} = nothing,
                      algorithm::JetAlgorithm.Algorithm = JetAlgorithm.AntiKt,
                      strategy = RecoStrategy.N2Tiled,
                      recombine = RecombinationMethods[RecombinationScheme.ESchemeRaw]) where {T <:
                                                                                               JetReconstruction.FourMomentum}
    Profile.init(n = 5 * 10^6, delay = 0.00001)
    function profile_events(events)
        for evt in events
            jet_reconstruct(evt; R = R, p = p, γ = γ, algorithm = algorithm,
                            strategy = strategy, recombine...)
        end
    end
    # Do a warm up run first to avoid JIT compilation costs
    profile_events(events)

    # Now take the actual profile
    @profile for i in 1:nsamples
        profile_events(events)
    end
    profile_dir = joinpath("profile", profile)
    mkpath(profile_dir)
    println("""Generating HTML flame graph at $(joinpath(profile_dir, "index.html"))""")
    statprofilehtml(path = profile_dir)
end

"""
    allocation_stats(events::Vector{Vector{T}}; algorithm::JetAlgorithm.Algorithm,
                          distance::Real = 0.4, p::Union{Real, Nothing} = nothing,
                          γ::Union{Real, Nothing} = nothing,
                          strategy::RecoStrategy.Strategy,
                          recombine = RecombinationMethods[RecombinationScheme.ESchemeRaw],
                          ptmin::Real = 5.0) where {T <: JetReconstruction.FourMomentum}
Take a memory allocation profile of the jet reconstruction code, printing the
output.
"""
function allocation_stats(events::Vector{Vector{T}}; algorithm::JetAlgorithm.Algorithm,
                          distance::Real = 0.4, p::Union{Real, Nothing} = nothing,
                          γ::Union{Real, Nothing} = nothing,
                          strategy::RecoStrategy.Strategy,
                          recombine = RecombinationMethods[RecombinationScheme.ESchemeRaw],
                          ptmin::Real = 5.0) where {T <: JetReconstruction.FourMomentum}
    println("Memory allocation statistics:")
    @timev for event in events
        _ = inclusive_jets(jet_reconstruct(event; R = distance, p = p, γ = γ,
                                           algorithm = algorithm,
                                           strategy = strategy, recombine...),
                           ptmin = ptmin)
    end
    nothing
end

"""
    benchmark_jet_reco(events::Vector{Vector{T}};
                            algorithm::JetAlgorithm.Algorithm,
                            distance::Real = 0.4,
                            p::Union{Real, Nothing} = nothing,
                            γ::Union{Real, Nothing} = nothing,
                            ptmin::Real = 5.0,
                            dcut = nothing,
                            njets = nothing,
                            strategy::RecoStrategy.Strategy,
                            nsamples::Integer = 1,
                            gcoff::Bool = false,
                            dump::Union{String, Nothing} = nothing,
                            dump_cs = false) where {T <: JetReconstruction.FourMomentum}

Benchmark the jet reconstruction code with given reconstruction parameters and
print summary statistics on the runtime.

# Notable Optional Arguments

- `dump`: If not `nothing`, write the list of reconstructed jets to a JSON
  formatted file.
- `dump_cs`: If `true`, dump the cluster sequence for each event.
"""
function benchmark_jet_reco(events::Vector{Vector{T}};
                            algorithm::JetAlgorithm.Algorithm,
                            distance::Real = 0.4,
                            p::Union{Real, Nothing} = nothing,
                            γ::Union{Real, Nothing} = nothing,
                            strategy::RecoStrategy.Strategy,
                            recombine = RecombinationMethods[RecombinationScheme.ESchemeRaw],
                            ptmin::Real = 5.0,
                            dcut = nothing,
                            njets = nothing,
                            nsamples::Integer = 1,
                            gcoff::Bool = false,
                            dump::Union{String, Nothing} = nothing,
                            dump_cs = false) where {T <: JetReconstruction.FourMomentum}

    # If we are dumping the results, setup the JSON structure
    if !isnothing(dump)
        jet_collection = FinalJets[]
    end

    # Vector for the trial results
    trial_timing = zeros(Float64, nsamples)

    # Set consistent algorithm power
    p = JetReconstruction.get_algorithm_power(p = p, algorithm = algorithm)
    @info "Jet reconstruction will use $(algorithm) with power $(p)"

    # Now setup timers and run the loop
    GC.gc()
    # Do a warm up run if we are running more than once
    start = nsamples > 1 ? 0 : 1
    for irun in start:nsamples
        gcoff && GC.enable(false)
        t_start = time_ns()
        for (ievt, event) in enumerate(events)
            cs = jet_reconstruct(event; R = distance, p = p, γ = γ, algorithm = algorithm,
                                 strategy = strategy, recombine...)
            if !isnothing(njets)
                selectedjets = exclusive_jets(cs; njets = njets)
            elseif !isnothing(dcut)
                selectedjets = exclusive_jets(cs; dcut = dcut)
            else
                selectedjets = inclusive_jets(cs; ptmin = ptmin)
            end
            # Only print the jet content once
            if irun == 0 || nsamples == 1
                @info begin
                    jet_output = "Event $(ievt)\n"
                    sort!(selectedjets, by = x -> pt(x), rev = true)
                    for (ijet, jet) in enumerate(selectedjets)
                        jet_output *= " $(ijet) - $(jet)\n"
                    end
                    "$(jet_output)"
                end
                if !isnothing(dump)
                    push!(jet_collection, FinalJets(ievt, final_jets(selectedjets)))
                end
                if dump_cs
                    println("Cluster sequence for event $(ievt)")
                    for (ijet, jet) in enumerate(cs.jets)
                        println(" $(ijet) - $(jet)")
                    end
                    for (ihistory, history) in enumerate(cs.history)
                        println(" $(ihistory) - $(history)")
                    end
                end
            end
        end
        t_stop = time_ns()
        gcoff && GC.enable(true)
        if irun > 0
            dt_μs = convert(Float64, t_stop - t_start) * 1.e-3
            trial_timing[irun] = dt_μs / length(events)
        end
    end

    if !isnothing(dump)
        open(dump, "w") do io
            JSON.print(io, jet_collection, 2)
        end
    end
    trial_timing
end

"""
    function print_statistics(trial_stats; outlier_exclusion = true, 
        outlier_band = 2, plot = false)

Give a statistical summary of trial runs. Statistics for the whole
sample are given, as well as those with outliers removed.

# Notes

## Why also record the lowest time? 
 
The argument is that on a "busy" machine, the run time of an application is
always TrueRunTime+Overheads, where Overheads is a nuisance parameter that
adds jitter, depending on the other things the machine is doing. Therefore
the minimum value is (a) more stable and (b) reflects better the intrinsic
code performance.
"""
function print_statistics(trial_stats; outlier_exclusion = true, outlier_band = 2.0,
                          plot = false)
    sstats = summarystats(trial_stats)
    iqr = sstats.q75 - sstats.q25
    println("Full statistics ($(sstats.nobs) samples)")
    println(pprint_trial_stats(sstats))

    if outlier_exclusion
        # The interquartile range (Q3 - Q1, or q75% - q25%) is the most useful way to filter outliers
        # as it's insensitive to the outliers themselves (unlike σ)
        min_val = sstats.q25 - outlier_band * iqr
        max_val = sstats.q75 + outlier_band * iqr
        no_outliers = trial_stats[(a -> a >= min_val && a <= max_val).(trial_stats)]
        no_outliers_stats = summarystats(no_outliers)
        println("Excluding outliers at $(outlier_band)xIQR (leaving $(no_outliers_stats.nobs) of $(sstats.nobs) samples)")
        println(pprint_trial_stats(no_outliers_stats))
    end

    if plot
        # Derive bin width from the IQR and number of samples
        bin_width = ceil(2iqr / sstats.nobs^(1 / 3))
        # No more than 80 bins (sometimes this gets way too big if there's an extreme outlier)
        hbins = min(ceil((sstats.max - sstats.min) / bin_width), 80)

        println(histogram(trial_stats, nbins = hbins, vertical = true,
                          title = "Histogram of event time per trial"))
        println(lineplot(collect(1:length(trial_stats)), trial_stats,
                         title = "Runtime per event across trials"))
    end
end

"""
    function pprint_trial_stats(sstats)

Generate a "pretty printed" string from a StatsBase statistics object,
which is basically use @sprintf to limit the precision
"""
function pprint_trial_stats(sstats)
    " - average time per event " * @sprintf("%.2f", sstats.mean) * " ± " *
    @sprintf("%.2f", sstats.sd) * " μs\n" *
    " - lowest time per event " * @sprintf("%.2f", sstats.min) * " μs"
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

        "--power", "--beta", "-p"
        help = """Power value for jet reconstruction (a.k.a. β for Valencia algorithm)"""
        arg_type = Float64

        "--gamma"
        help = """Gamma, γ, only for Valencia algorithm)"""
        arg_type = Float64

        "--strategy", "-S"
        help = """Strategy for the algorithm, valid values: $(join(JetReconstruction.AllJetRecoStrategies, ", "))"""
        arg_type = RecoStrategy.Strategy
        default = RecoStrategy.Best

        "--recombine"
        help = """Recombination scheme to use for jet reconstruction: $(join(JetReconstruction.AllRecombinationSchemes, ", "))"""
        arg_type = RecombinationScheme.Recombine
        default = RecombinationScheme.ESchemeRaw

        "--nsamples", "-m"
        help = "Number of measurement points to acquire."
        arg_type = Int
        default = 1

        "--dump-clusterseq"
        help = "Dump the cluster sequence for each event"
        action = :store_true

        "--gcoff"
        help = "Turn off Julia garbage collector during each time measurement."
        action = :store_true

        "--profile"
        help = """Profile code and generate a flame graph; the results 
        will be stored in 'profile/PROFILE' subdirectory, where PROFILE is the value
        given to this option)"""

        "--alloc"
        help = "Provide memory allocation statistics."
        action = :store_true

        "--plot"
        help = "Plot a histogram of trail times on the terminal"
        action = :store_true

        "--dump"
        help = "Write list of reconstructed jets to a JSON formatted file"

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

function main()
    args = parse_command_line(ARGS)
    if args[:debug]
        logger = ConsoleLogger(stdout, Logging.Debug)
    elseif args[:info]
        logger = ConsoleLogger(stdout, Logging.Info)
    else
        logger = ConsoleLogger(stdout, Logging.Warn)
    end
    global_logger(logger)

    if isnothing(args[:algorithm]) && isnothing(args[:power])
        @warn "Neither algorithm nor power specified, defaulting to AntiKt"
        args[:algorithm] = JetAlgorithm.AntiKt
    end

    # Try to read events into the correct type!
    if JetReconstruction.is_ee(args[:algorithm])
        jet_type = EEJet
    else
        jet_type = PseudoJet
    end
    events::Vector{Vector{jet_type}} = read_final_state_particles(args[:file], jet_type;
                                                                  maxevents = args[:maxevents],
                                                                  skipevents = args[:skip])

    # Major switch between modes of running
    trial_stats = nothing
    if args[:alloc]
        allocation_stats(events; distance = args[:distance],
                         p = args[:power], γ = args[:gamma], algorithm = args[:algorithm],
                         strategy = args[:strategy],
                         recombine = JetReconstruction.RecombinationMethods[args[:recombine]],
                         ptmin = args[:ptmin])
    elseif !isnothing(args[:profile])
        profile_code(events, args[:profile], args[:nsamples];
                     R = args[:distance], p = args[:power], γ = args[:gamma],
                     algorithm = args[:algorithm], strategy = args[:strategy],
                     recombine = JetReconstruction.RecombinationMethods[args[:recombine]])
    else
        trial_stats = benchmark_jet_reco(events, distance = args[:distance],
                                         algorithm = args[:algorithm],
                                         p = args[:power], γ = args[:gamma],
                                         strategy = args[:strategy],
                                         recombine = JetReconstruction.RecombinationMethods[args[:recombine]],
                                         ptmin = args[:ptmin], dcut = args[:exclusive_dcut],
                                         njets = args[:exclusive_njets],
                                         nsamples = args[:nsamples], gcoff = args[:gcoff],
                                         dump = args[:dump],
                                         dump_cs = args[:dump_clusterseq])
    end
    if !isnothing(trial_stats)
        print_statistics(trial_stats; plot = args[:plot])
    end
    nothing
end

main()
