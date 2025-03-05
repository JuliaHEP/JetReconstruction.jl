
using ArgParse
using Profile
using Logging
using JSON

using LorentzVectorHEP
using JetReconstruction
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

function set_up()
#main 
end 

#
function select_ABS_RAP_max(event, absrapmax)
    filter(event -> abs(event.rapidity) <= absrapmax, event)
end

function main()
    print("Zdravei")

    hard_event = PseudoJet[]
    full_event = PseudoJet[]    

    rapmax = 5.0 
    hard_event = select_ABS_RAP_max(hard_event,rapmax)
    full_event = select_ABS_RAP_max(full_event,rapmax)

    args = parse_command_line(ARGS)
    logger = ConsoleLogger(stdout, Logging.Info)
    global_logger(logger)

    #if JetReconstruction.is_ee(args[:algorithm])
    #    jet_type = EEjet
    #else
       jet_type = PseudoJet
    #end

    events::Vector{Vector{jet_type}} = read_final_state_particles(args[:file],
                                        maxevents = args[:maxevents],
                                        skipevents = args[:skip],
                                        T = jet_type)
    #if isnothing(args[:algorithm]) && isnothing(args[:power])
    #    @warn "Neither algorithm nor power specified, defaulting to AntiKt"
    #    args[:algorithm] = JetAlgorithm.AntiKt
    #end    

    grid_size = 0.4
    soft_killer_grid = SoftKiller(rapmax, grid_size)

    #helper that does everthying from main in example.cc
                                       
end

main()