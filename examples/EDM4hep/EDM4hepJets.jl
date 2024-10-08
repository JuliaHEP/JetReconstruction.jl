using ArgParse
using EDM4hep
using EDM4hep.RootIO
using JetReconstruction
using Logging

function parse_command_line(args)
    s = ArgParseSettings(autofix_names = true)
    @add_arg_table! s begin
        "--maxevents", "-n"
        help = "Maximum number of events to read (default 1; -1 means read all events)."
        arg_type = Int
        default = 1

        "file"
        help = "EDM4hep events file to read"
        required = true
    end
    return parse_args(args, s; as_symbols = true)
end

function main()
    args = parse_command_line(ARGS)
    logger = ConsoleLogger(stdout, Logging.Info)
    global_logger(logger)

    # Open input event file
    reader = RootIO.Reader(args[:file])
    events = RootIO.get(reader, "events")

    # Reconstruct each event
    for (ievt, evt) in enumerate(events)
        if args[:maxevents] > 0 && ievt > args[:maxevents]
            break
        end
        recps = RootIO.get(reader, evt, "ReconstructedParticles")
        particles = ReconstructedParticle[]
        for recp in recps
            push!(particles, recp)
        end

        cs = jet_reconstruct(particles; algorithm = JetAlgorithm.Durham)
        @info begin
            jets = "Event $(ievt)\n"
            for jet in exclusive_jets(cs; njets = 2, T = EEjet)
                jets *= " $jet\n"
            end
            jets
        end
    end
end

main()
