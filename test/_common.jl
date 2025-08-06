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
using StructArrays

using JetReconstruction: EERecoJet,
                         fill_reco_array!, insert_new_jet!, copy_to_slot!,
                         dij_dist, valencia_distance

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

"""
    struct ComparisonTestValencia

Test parameters for Valencia jet algorithm comparison, with explicit gamma parameter.
"""
struct ComparisonTestValencia
    events_file::AbstractString
    fastjet_outputs::AbstractString
    algorithm::JetAlgorithm.Algorithm
    strategy::RecoStrategy.Strategy
    power::Real
    R::Real
    γ::Real
    selector::Any
    selector_name::AbstractString
    recombine::Any
    preprocess::Any
end

function ComparisonTestValencia(events_file, fastjet_outputs, algorithm, strategy, power, R,
                                γ, selector, selector_name)
    ComparisonTestValencia(events_file, fastjet_outputs, algorithm, strategy, power, R, γ,
                           selector, selector_name, addjets, nothing)
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

function run_reco_test(test::ComparisonTest; testname = nothing,
                       break_history_indices = false)
    jet_type = JetReconstruction.is_ee(test.algorithm) ? EEJet : PseudoJet
    # Read the input events
    events = JetReconstruction.read_final_state_particles(test.events_file, jet_type)
    # mess around with indices if requested
    if break_history_indices
        for event in events
            event = [jet_type(jet; cluster_hist_index = 0) for jet in event]
        end
    end
    # Read the fastjet results
    fastjet_jets = read_fastjet_outputs(test.fastjet_outputs)
    sort_jets!(fastjet_jets)

    # Run the jet reconstruction
    jet_collection = Vector{FinalJets}()
    for (ievent, event) in enumerate(events)
        if test.algorithm == JetAlgorithm.Valencia
            # For VLC: pass both beta (power) and γ
            cluster_seq = JetReconstruction.jet_reconstruct(event; R = test.R,
                                                            p = test.power,
                                                            γ = getfield(test, :γ),
                                                            algorithm = test.algorithm,
                                                            strategy = test.strategy,
                                                            recombine = test.recombine,
                                                            preprocess = test.preprocess)
        else
            # All other algorithms
            cluster_seq = JetReconstruction.jet_reconstruct(event; R = test.R,
                                                            p = test.power,
                                                            algorithm = test.algorithm,
                                                            strategy = test.strategy,
                                                            recombine = test.recombine,
                                                            preprocess = test.preprocess)
        end
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
                        @test jet.pt≈fastjet_jets[ievt]["jets"][ijet]["pt"] rtol=1e-7
                    end
                end
            end
        end
    end
end

function run_reco_test(test::ComparisonTestValencia; testname = nothing)
    events = JetReconstruction.read_final_state_particles(test.events_file)
    fastjet_jets = read_fastjet_outputs(test.fastjet_outputs)
    sort_jets!(fastjet_jets)

    jet_collection = Vector{FinalJets}()
    for (ievent, event) in enumerate(events)
        cluster_seq = JetReconstruction.jet_reconstruct(event; R = test.R, p = test.power,
                                                        γ = test.γ,
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
        for (ievt, event) in enumerate(jet_collection)
            @testset "Event $(ievt)" begin
                @test size(event.jets) == size(fastjet_jets[ievt]["jets"])
                # For each reconstructed jet, find the closest reference jet in deltaR
                ref_jets = fastjet_jets[ievt]["jets"]
                used_refs = falses(length(ref_jets))
                function deltaR(rap1, phi1, rap2, phi2)
                    dphi = abs(phi1 - phi2)
                    if dphi > π
                        dphi -= 2π
                    end
                    return sqrt((rap1 - rap2)^2 + dphi^2)
                end
                for (ijet, jet) in enumerate(event.jets)
                    # Find closest reference jet in deltaR that hasn't been used
                    minidx = 0
                    mindr = Inf
                    norm_phi = jet.phi < 0.0 ? jet.phi + 2π : jet.phi
                    for (ridx, rjet) in enumerate(ref_jets)
                        if used_refs[ridx]
                            continue
                        end
                        norm_phi_ref = rjet["phi"] < 0.0 ? rjet["phi"] + 2π : rjet["phi"]
                        dr = deltaR(jet.rap, norm_phi, rjet["rap"], norm_phi_ref)
                        if dr < mindr
                            mindr = dr
                            minidx = ridx
                        end
                    end
                    if minidx == 0
                        println("No unused reference jet found for reconstructed jet $(ijet)")
                        continue
                    end
                    used_refs[minidx] = true
                    rjet = ref_jets[minidx]
                    rap_ref = rjet["rap"]
                    phi_ref = rjet["phi"]
                    pt_ref = rjet["pt"]
                    normalised_phi_ref = phi_ref < 0.0 ? phi_ref + 2π : phi_ref
                    rap_test = isapprox(jet.rap, rap_ref; atol = 1e-4)
                    phi_test = isapprox(norm_phi, normalised_phi_ref; atol = 1e-4)
                    pt_test = isapprox(jet.pt, pt_ref; rtol = 1e-5)
                    if !rap_test || !phi_test || !pt_test
                        println("Jet mismatch in Event $(ievt), Jet $(ijet):")
                        println("  Failing Jet: pt=$(jet.pt), rap=$(jet.rap), phi=$(norm_phi)")
                        println("  Reference Jet: pt=$(pt_ref), rap=$(rap_ref), phi=$(phi_ref)")
                        println("  Passes: pt=$(pt_test), rap=$(rap_test), phi=$(phi_test)")
                        println("\nAll jets in this event:")
                        for (jidx, j) in enumerate(event.jets)
                            norm_phi_j = j.phi < 0.0 ? j.phi + 2π : j.phi
                            println("  Jet $(jidx): pt=$(j.pt), rap=$(j.rap), phi=$(norm_phi_j)")
                        end
                        println("\nAll reference jets in this event:")
                        for (ridx, rjet2) in enumerate(ref_jets)
                            println("  Ref Jet $(ridx): pt=$(rjet2["pt"]), rap=$(rjet2["rap"]), phi=$(rjet2["phi"])")
                        end
                    end
                    @test rap_test
                    @test phi_test
                    @test pt_test
                end
            end
        end
    end
end
