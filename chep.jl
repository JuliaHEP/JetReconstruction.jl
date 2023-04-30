#! /usr/bin/env julia
"""
Wrapper to run Atell's jet reco code feeding in the standard set of HepMC events that
are used for benchmarking jet reco.
"""

import HepMC3
using ArgParse

using JetReconstruction

const R = 0.4
const ptmin = 5.0

read_events(fname, maxevents=-1, skipevents = 0) = begin
    f = open(fname)

    events = Vector{PseudoJet}[]

    ipart = 1
    HepMC3.read_events(f, maxevents=maxevents, skipevents=skipevents) do parts
        input_particles = PseudoJet[]
        for p in parts
            if p.status == 1
                push!(input_particles, PseudoJet(p.momentum.px, p.momentum.py, p.momentum.pz, p.momentum.e))
            end
        end
        push!(events, input_particles)
        println("PseudoJet number $ipart: ", length(input_particles))
        ipart += 1
    end

    println("Total Events: ", length(events))
    events
end

final_jets(jets::Vector{Vector{Float64}}, ptmin::AbstractFloat) = begin
    count = 0
    for jet in jets
        dcut = ptmin*ptmin
        p = PseudoJet(jet[1], jet[2], jet[3], jet[4])
        if p._pt2 > dcut
            count += 1
            println(p)
            # println("$(count), $(rap(p)), $(phi(p)), $(pt2(p))")
        end
    end
end
 
in_mem_process(events::Vector{Vector{PseudoJet}}) = begin
    println("Will process $(size(events)) events")
    for event in events
        # Perhaps inefficient to use Vector of Vectors? (Matrix instead?)
        jet_struct = Vector{Vector{Float64}}(undef, size(event)[1])
        for (ipart, initial_particle) in enumerate(event)
            particle_data = [initial_particle.px, initial_particle.py, initial_particle.pz, initial_particle.E]
            jet_struct[ipart] = particle_data
        end
        smalljets, smallind = anti_kt_algo(jet_struct, R=0.4)
        # println(smalljets, smallind)
        println(typeof(smalljets))
        final_jets(smalljets, ptmin)
    end
end

parse_command_line(args) = begin
    s = ArgParseSettings(autofix_names=true)
    @add_arg_table s begin
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
        "--profile-only"
        help = "Profile code and generate a flame graph. Skip the code timing."
        action = :store_true
        "--alloc"
        help = "Provide memory allocation statistics."
        action = :store_true
        "--dump"
        help = "Display list of recontructed jets"
        action = :store_true
        "--json"
        help = "Write list of recontructed jets to a JSON file"
        action = :store_true
        "file"
        help = "HepMC3 event file in HepMC3 to read."
        required = true
    end
    return parse_args(args, s; as_symbols = true)
end

main() = begin
    args = parse_command_line(ARGS)
    events::Vector{Vector{PseudoJet}} = read_events(args[:file], args[:maxevents], args[:skip])
    in_mem_process(events)
    nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
