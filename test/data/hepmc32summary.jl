#! /usr/bin/env julia

using ArgParse
using JetReconstruction
using LorentzVectorHEP
using Statistics
using UnicodePlots

function parse_command_line(args)
    s = ArgParseSettings(autofix_names = true)
    @add_arg_table! s begin
        "--summary"
        help = "Print only summary information, filename and average density"
        action = :store_true

        "--maxevents", "-n"
        help = "Maximum number of events to read. -1 to read all events from the  file."
        arg_type = Int
        default = -1

        "--skip", "-s"
        help = "Number of events to skip at beginning of the file."
        arg_type = Int
        default = 0

        "files"
        help = "The HepMC3 event files to read."
        required = true
        nargs = '+'
    end
    return parse_args(args, s; as_symbols = true)
end

function main()
    args = parse_command_line(ARGS)

    for file in args[:files]
        events = read_final_state_particles(file, LorentzVector{Float64};
                                            maxevents = args[:maxevents],
                                            skipevents = args[:skip])
        n_events = length(events)
        n_particles = Int[]
        for e in events
            push!(n_particles, length(e))
        end
        average_n = mean(n_particles)
        if args[:summary]
            println("$file,$average_n")
        else
            println("File $file")
            println("  Number of events: $n_events")
            println("  Average number of particles: ", mean(n_particles))
            if n_events > 1
                println(histogram(n_particles))
            end
        end
    end
end

main()
