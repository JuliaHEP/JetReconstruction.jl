#! /usr/bin/env julia

"""
Process independent events concurrently with one reusable N2Tiled workspace
owned by each long-lived worker task.

Run with, for example:

```sh
julia --threads=auto --project=examples examples/n2tiled-multithreaded.jl \
    --maxevents=100 test/data/events.pp13TeV.hepmc3.zst
```
"""

using ArgParse
using JetReconstruction
using LorentzVectorHEP

"""
    threaded_n2tiled_reconstruction(events; R=0.4, ptmin=5.0)

Reconstruct `events` using dynamically scheduled worker tasks. Each worker owns
one `N2TiledWorkspace`; only independently owned inclusive-jet selections are
stored after the callback returns.
"""
function threaded_n2tiled_reconstruction(events;
                                         R::Real = 0.4,
                                         ptmin::Real = 5.0)
    selected_jets = Vector{Vector{LorentzVector{Float64}}}(undef, length(events))
    isempty(events) && return selected_jets

    nworkers = min(Threads.nthreads(), length(events))

    # A shared channel balances events dynamically. Workspaces belong to the
    # tasks below rather than to thread IDs, so task migration is safe.
    jobs = Channel{Int}(length(events))
    for event_index in eachindex(events)
        put!(jobs, event_index)
    end
    close(jobs)

    @sync for _ in 1:nworkers
        Threads.@spawn begin
            workspace = N2TiledWorkspace()

            for event_index in jobs
                selected_jets[event_index] = with_n2tiled_reconstruction(workspace,
                                                                         events[event_index];
                                                                         algorithm = JetAlgorithm.AntiKt,
                                                                         R = R,) do clusterseq
                    # `inclusive_jets` returns an independently owned vector.
                    # Use `deepcopy(clusterseq)` here instead when the complete
                    # clustering sequence must outlive this callback.
                    inclusive_jets(clusterseq; ptmin = ptmin)
                end
            end
        end
    end

    return selected_jets
end

function parse_command_line(args)
    settings = ArgParseSettings(autofix_names = true)
    @add_arg_table! settings begin
        "--maxevents", "-n"
        help = "Maximum number of events to read; -1 reads all events."
        arg_type = Int
        default = -1

        "--ptmin"
        help = "Minimum transverse momentum for inclusive jets."
        arg_type = Float64
        default = 5.0

        "--distance", "-R"
        help = "Jet radius parameter."
        arg_type = Float64
        default = 0.4

        "file"
        help = "HepMC3 event file to read."
        required = true
    end

    return parse_args(args, settings; as_symbols = true)
end

function main(args = ARGS)
    options = parse_command_line(args)
    events = read_final_state_particles(options[:file], PseudoJet;
                                        maxevents = options[:maxevents])

    selected_jets = threaded_n2tiled_reconstruction(events;
                                                    R = options[:distance],
                                                    ptmin = options[:ptmin])

    println("Processed $(length(events)) events with $(Threads.nthreads()) Julia threads; " *
            "selected $(sum(length, selected_jets)) jets.")

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
