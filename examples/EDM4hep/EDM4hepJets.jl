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

        "--printjets"
        help = "Print reconstructed jets (selecting njets from the reconstructed event)"
        action = :store_true

        "--njets"
        help = "Number of jets to output"
        arg_type = Int
        default = 2

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
    global_timer_start = time_ns()
    reader = RootIO.Reader(args[:file])
    events = RootIO.get(reader, "events")

    total_time = 0.0
    count = 0

    # Reconstruct each event
    for (ievt, evt) in enumerate(events)
        count += 1
        if args[:maxevents] > 0 && ievt > args[:maxevents]
            break
        end
        recps = RootIO.get(reader, evt, "ReconstructedParticles")

        start_jet = time_ns()
        cs = jet_reconstruct(recps; algorithm = JetAlgorithm.Durham)
        end_jet = time_ns()
        total_time += end_jet - start_jet
        if args[:printjets]
            @info begin
                jets = "Event $(ievt)\n"
                for jet in exclusive_jets(cs, EEJet; njets = args[:njets])
                    jets *= " $jet\n"
                end
                jets
            end
        end
    end
    global_timer_stop = time_ns()
    @info "Reconstructed $(count) events in $(total_time/1.0e9)s: rate = $((count * 1.0e9) / total_time)Hz"
    @info "Total time: $((global_timer_stop - global_timer_start)/1.0e9)s"
    count, total_time
end

main()
