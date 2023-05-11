#! /usr/bin/env julia
"""
Wrapper to run jet reco code feeding in the standard set of HepMC events that
are used for benchmarking jet reco.
"""

import HepMC3
import FlameGraphs
import ProfileSVG

using ArgParse
using Profile
using Colors
using StatProfilerHTML

using JetReconstruction

const R = 0.4
const ptmin = 5.0

function read_events(fname; maxevents=-1, skipevents=0)
    f = open(fname)

    events = Vector{PseudoJet}[]

    ipart = 1
    HepMC3.read_events(f, maxevents=maxevents, skipevents=skipevents) do parts
        input_particles = PseudoJet[]
        for p in parts
            if p.status == 1
                push!(
                    input_particles,
                    PseudoJet(p.momentum.px, p.momentum.py, p.momentum.pz, p.momentum.e),
                )
            end
        end
        push!(events, input_particles)
        ipart += 1
    end

    println("Total Events: ", length(events))
    events
end

function pseudojets2vectors(events::Vector{Vector{PseudoJet}})
    """Convert events in PseudoJets into deep vectors"""
    event_vector = Vector{Vector{Vector{Float64}}}(undef, size(events)[1])
    for (ievent, event) in enumerate(events)
        jet_struct = Vector{Vector{Float64}}(undef, size(event)[1])
        for (ipart, initial_particle) in enumerate(event)
            jet_struct[ipart] = [
                initial_particle.px,
                initial_particle.py,
                initial_particle.pz,
                initial_particle.E,
            ]
        end
        event_vector[ievent] = jet_struct
    end
    event_vector
end

function final_jets(jets::Vector{Vector{Float64}}, ptmin::AbstractFloat)
    count = 0
    final_jets = Vector{FinalJet}()
    sizehint!(final_jets, 6)
    for jet in jets
        dcut = ptmin^2
        p = PseudoJet(jet[1], jet[2], jet[3], jet[4])
        if p._pt2 > dcut
            count += 1
            push!(final_jets, FinalJet(rap(p), phi(p), sqrt(pt2(p))))
        end
    end
    final_jets
end

function profile_code(events, niters)
    Profile.init(n=10^6, delay=0.00001)
    profile_events(events) = begin
        for evt in events
            anti_kt_algo(evt, R=0.4)
        end
    end
    profile_events(events[1:1])
    @profile for i = 1:niters
        profile_events(events)
    end
    statprofilehtml()
    fcolor = FlameGraphs.FlameColors(
        reverse(colormap("Blues", 15))[1:5],
        colorant"slategray4",
        colorant"gray95",
        reverse(colormap("Reds", 15))[1:5],
        reverse(sequential_palette(39, 10; s=38, b=2))[1:5],#brownish pallette
    )
    ProfileSVG.save(
        fcolor,
        joinpath("statprof", "profsvg.svg");
        combine=true,
        timeunit=:ms,
        font="Arial, Helvetica, sans-serif"
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
    nsamples::Integer=1,
    gcoff::Bool=false,
    profile::Bool=false,
    dump::Union{String,Nothing}=nothing
)
    println("Will process $(size(events)) events")

    # First, convert all events into the Vector of Vectors that Atell's
    # code likes
    event_vector = pseudojets2vectors(events)

    # If we are dumping the results, setup the JSON structure
    if !isnothing(dump)
        jet_collection = FinalJets[]
    end

    # Warmup code if we are doing a multi-sample timing run
    if nsamples > 1 || profile
        println("Doing initial warm-up run")
        for event in event_vector
            anti_kt_algo(event, R=0.4)
        end
    end

    GC.gc()
    gcoff && GC.enable(false)

    if profile
        profile_code(event_vector, nsamples)
        return Nothing
    end

    # Now setup timers and run the loop
    cummulative_time = 0.0
    cummulative_time2 = 0.0
    for irun = 1:nsamples
        print("$(irun)/$(nsamples) ")
        t_start = time_ns()
        for (ievt, event) in enumerate(event_vector)
            finaljets, _ = anti_kt_algo(event, R=0.4)
            fj = final_jets(finaljets, ptmin)
            if !isnothing(dump) && irun == 1
                println("Event $(ievt)")
                for (ijet, jet) in enumerate(fj)
                    println(" $(ijet) - $(jet)")
                end
                push!(jet_collection, FinalJets(ievt, fj))
            end
        end
        t_stop = time_ns()
        dt_μs = convert(Float64, t_stop - t_start) * 1.e-3
        println(dt_μs)
        cummulative_time += dt_μs
        cummulative_time2 += dt_μs^2
    end

    gcoff && GC.enable(true)

    mean = cummulative_time / nsamples
    cummulative_time2 /= nsamples
    if nsamples > 1
        sigma = sqrt(nsamples / (nsamples - 1) * (cummulative_time2 - mean^2))
    else
        sigma = 0.0
    end
    mean /= length(events)
    sigma /= length(events)
    println("Processed $(length(events)) events $(nsamples) times")
    println("Time per event $(mean) ± $(sigma) μs")

    if !isnothing(dump)
        open(dump, "w") do io
            JSON3.pretty(io, jet_collection)
        end
    end
end

parse_command_line(args) = begin
    s = ArgParseSettings(autofix_names=true)
    @add_arg_table! s begin
        "--maxevents", "-n"
        help = "Maximum number of events to read. -1 to read all events from the  file."
        arg_type = Int
        default = -1

        "--skip", "-s"
        help = "Number of events to skip at beginning of the file."
        arg_type = Int
        default = 0

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

        "--dump"
        help = "Write list of recontructed jets to a JSON formatted file and also to stdout"

        "file"
        help = "HepMC3 event file in HepMC3 to read."
        required = true
    end
    return parse_args(args, s; as_symbols=true)
end

main() = begin
    args = parse_command_line(ARGS)
    events::Vector{Vector{PseudoJet}} =
        read_events(args[:file], maxevents=args[:maxevents], skipevents=args[:skip])
    jet_process(events, nsamples=args[:nsamples], gcoff=args[:gcoff], profile=args[:profile],
        dump=args[:dump])
    nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
