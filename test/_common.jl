# Put here common parameters for the jet reconstruction tests
# This allows different test scripts to be run standalone, but
# to easily share common parameters.

# Set this variable for use by the include guard wrapper
const JETRECO_TEST_COMMON = true

using JetReconstruction
using Logging
using LorentzVectorHEP
using JSON
using Test

logger = ConsoleLogger(stdout, Logging.Warn)
global_logger(logger)

const events_file_pp = joinpath(@__DIR__, "data", "events.pp13TeV.hepmc3.zst")
const events_file_ee = joinpath(@__DIR__, "data", "events.eeH.hepmc3.zst")

const pp_algorithms = Dict(-1 => "Anti-kt",
                           0 => "Cambridge/Aachen",
                           1 => "Inclusive-kt",
                           1.5 => "Generalised-kt")

"""
    struct ComparisonTest

Test parameters for the comparison of jet reconstruction against known FastJet
results.

# Fields
- `events_file::AbstractString`: The file containing the event data.
- `fastjet_outputs::AbstractString`: The file containing the FastJet outputs in
  JSON format (N.B. can be gzipped)
- `algorithm::JetAlgorithm.Algorithm`: The algorithm used for jet
  reconstruction.
- `power::Real`: The power parameter for the jet reconstruction algorithm.
- `R::Real`: The radius parameter for the jet reconstruction algorithm.
- `selector`: The selector used for final jet outputs.

`selector` should be a closure that takes a `ClusterSequence` object and returns
a vector of `FinalJet` objects, e.g.,
```julia
selector = (cs) -> exclusive_jets(cs; njets = 2)
```
"""
struct ComparisonTest
    events_file::AbstractString
    fastjet_outputs::AbstractString
    algorithm::JetAlgorithm.Algorithm
    strategy::RecoStrategy.Strategy
    power::Real
    R::Real
    selector::Any
    selector_name::AbstractString
    recombine::Any
    preprocess::Any
end

"""Constructor where there is no selector_name given"""
function ComparisonTest(events_file, fastjet_outputs, algorithm, strategy, power, R,
                        selector)
    selector_name = ""
    ComparisonTest(events_file, fastjet_outputs, algorithm, strategy, power, R, selector,
                   selector_name, addjets, nothing)
end

"""Constructor with no recombine or preprocess specified"""
function ComparisonTest(events_file, fastjet_outputs, algorithm, strategy, power, R,
                        selector, selector_name)
    ComparisonTest(events_file, fastjet_outputs, algorithm, strategy, power, R, selector,
                   selector_name, addjets, nothing)
end

"""Read JSON file with fastjet jets in it"""
function read_fastjet_outputs(fname)
    f = JetReconstruction.open_with_stream(fname)
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

function sort_jets!(jet_array::Vector{<:LorentzVectorCyl})
    jet_pt(jet) = jet.pt
    sort!(jet_array, by = jet_pt, rev = true)
end

function run_reco_test(test::ComparisonTest; testname = nothing)
    # Read the input events
    events = JetReconstruction.read_final_state_particles(test.events_file)
    # Read the fastjet results
    fastjet_jets = read_fastjet_outputs(test.fastjet_outputs)
    sort_jets!(fastjet_jets)

    # Run the jet reconstruction
    jet_collection = Vector{FinalJets}()
    for (ievent, event) in enumerate(events)
        cluster_seq = JetReconstruction.jet_reconstruct(event; R = test.R, p = test.power,
                                                        algorithm = test.algorithm,
                                                        strategy = test.strategy,
                                                        recombine = test.recombine,
                                                        preprocess = test.preprocess)
        finaljets = final_jets(test.selector(cluster_seq))
        sort_jets!(finaljets)
        push!(jet_collection, FinalJets(ievent, finaljets))
    end

    if isnothing(testname)
        testname = "FastJet comparison: alg=$(test.algorithm), p=$(test.power), R=$(test.R), strategy=$(test.strategy)"
        if test.selector_name != ""
            testname *= ", $(test.selector_name)"
        end
    end

    @testset "$testname" begin
        # Test each event in turn...
        for (ievt, event) in enumerate(jet_collection)
            @testset "Event $(ievt)" begin
                @test size(event.jets) == size(fastjet_jets[ievt]["jets"])
                # Test each jet in turn
                for (ijet, jet) in enumerate(event.jets)
                    if ijet <= size(fastjet_jets[ievt]["jets"])[1]
                        # Approximate test - note that @test macro passes the 
                        # tolerance into the isapprox() function
                        # Use atol for position as this is absolute, but rtol for
                        # the momentum
                        # Sometimes phi could be in the range [-π, π], but FastJet always is [0, 2π]
                        normalised_phi = jet.phi < 0.0 ? jet.phi + 2π : jet.phi
                        @test jet.rap≈fastjet_jets[ievt]["jets"][ijet]["rap"] atol=1e-7
                        @test normalised_phi≈fastjet_jets[ievt]["jets"][ijet]["phi"] atol=1e-7
                        @test jet.pt≈fastjet_jets[ievt]["jets"][ijet]["pt"] rtol=1e-6
                    end
                end
            end
        end
    end
end
