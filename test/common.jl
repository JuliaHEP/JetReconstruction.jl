# Put here common parameters for the jet reconstruction tests
# This allows different test scripts to be run standalone, but
# to easily share common parameters.

const events_file_pp = joinpath(@__DIR__, "data", "events.pp13TeV.hepmc3.gz")
const events_file_ee = joinpath(@__DIR__, "data", "events.eeH.hepmc3.gz")

const pp_algorithms = Dict(-1 => "Anti-kt",
                        0 => "Cambridge/Achen",
                        1 => "Inclusive-kt",
                        1.5 => "Generalised-kt")

"""Simple structure with necessary parameters for an exclusive selection test"""
struct InclusiveTest
    algname::AbstractString
    selction::AbstractString
    power::Int
    fastjet_file::AbstractString
    dijmax::Any
    njets::Any
    function InclusiveTest(selction, power, fastjet_file, dijmax, njets)
        new(pp_algorithms[power], selction, power, fastjet_file, dijmax, njets)
    end
end

"""Read JSON file with fastjet jets in it"""
function read_fastjet_outputs(fname)
    f = open(fname)
    JSON.parse(f)
end

"""Sort jet outputs by pt of final jets"""
function sort_jets!(event_jet_array)
    jet_pt(jet) = jet["pt"]
    for e in event_jet_array
        sort!(e["jets"], by = jet_pt, rev = true)
    end
end

function sort_jets!(jet_array::Vector{FinalJet})
    jet_pt(jet) = jet.pt
    sort!(jet_array, by = jet_pt, rev = true)
end

function sort_jets!(jet_array::Vector{LorentzVectorCyl})
    jet_pt(jet) = jet.pt
    sort!(jet_array, by = jet_pt, rev = true)
end