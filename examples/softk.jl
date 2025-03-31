
using ArgParse
using Profile
using Logging
using JSON

using LorentzVectorHEP
using JetReconstruction
using Plots 
using Pkg
Pkg.develop(path="/Users/emadimtrova/Desktop/uni/spring25/UTRA/JetReconstruction.jl")

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

        "--algorithm", "-A"
        help = """Algorithm to use for jet reconstruction: $(join(JetReconstruction.AllJetRecoAlgorithms, ", "))"""
        arg_type = JetAlgorithm.Algorithm

        "--power", "-p"
        help = """Power value for jet reconstruction"""
        arg_type = Float64
        
        #add arguments if needed

        "file"
        help = "HepMC3 event file in HepMC3 to read."
        required = true
    end
    return parse_args(args, s; as_symbols = true)
end

read_event(fname; maxevents = -1, skipevents = 0, T = PseudoJet ) = begin
    hard_event = PseudoJet[]
    full_event = PseudoJet[]    
    events = read_final_state_particles(fname; maxevents, skipevents, T)

    #plot_graph(full_event, hard_event,events)
    Y = Float64[]
    Phi = Float64[]
    pt = Float64[]
    nsub = 0 

    for ev in events
        input_p = Vector{PseudoJet}()
        push_data(ev, input_p, Y, Phi,pt);
        append!(full_event, input_p) 
        nsub += 1

    end

    if (nsub == 1) #in the case of only one event after the reading 
        hard_event .= full_event 
        nsub += 1
    end
    if (nsub == 0)
        throw("Error: read empty event\n")
    end 

    plot_set_up(Y, Phi, pt, "Event before SoftKiller")

    return hard_event, full_event
end

function main()
    
    args = parse_command_line(ARGS)
    logger = ConsoleLogger(stdout, Logging.Info)
    global_logger(logger)

    #if JetReconstruction.is_ee(args[:algorithm])
    #    jet_type = EEjet
    #else
       jet_type = PseudoJet
    #end

    hard_event = PseudoJet[]
    full_event = PseudoJet[]   
    reduced_event = PseudoJet[] 

    #read event 

    hard_event, full_event = read_event(args[:file],
                                        maxevents = args[:maxevents],
                                        skipevents = args[:skip],
                                        T = jet_type)

    rapmax = 5.0 

    hard_event = select_ABS_RAP_max(hard_event,rapmax)  
    full_event = select_ABS_RAP_max(full_event,rapmax)
    
    #clustering

    grid_size = 0.4
    soft_killer = SoftKiller(rapmax, grid_size)

    pt_threshold = 0.00
    soft_killer_event = PseudoJet[]   

    reduced_event, pt_threshold = apply(soft_killer, full_event, soft_killer_event, pt_threshold)
    print("Pth value: ", pt_threshold, "\n")
end

main()