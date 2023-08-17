#! /usr/bin/env julia

using JetReconstruction
using Test
using JSON

"""Read JSON file with fastjet jets in it"""
function read_fastjet_outputs(;fname="test/data/jet_collections_fastjet.json")
    f = open(fname)
    JSON.parse(f)
end

"""Sort jet outputs by pt of final jets"""
function sort_jets!(event_jet_array)
    jet_pt(jet) = jet["pt"]
    for e in event_jet_array
        sort!(e["jets"], by=jet_pt, rev=true)
    end
end

function sort_jets!(jet_array::Vector{FinalJet})
    jet_pt(jet) = jet.pt
    sort!(jet_array, by=jet_pt, rev=true)
end

function main()
    # Read our fastjet outputs
    fastjet_jets = read_fastjet_outputs()
    sort_jets!(fastjet_jets)

    # Test each stratgy...
    do_jet_test(N2Plain, fastjet_jets)
    do_jet_test(N2TiledSoAGlobal, fastjet_jets)
    do_jet_test(N2TiledSoATile, fastjet_jets)
    do_jet_test(N2TiledLL, fastjet_jets)

    # Atell's original test
    original_tests()
end

function do_jet_test(strategy::JetRecoStrategy, fastjet_jets;
    ptmin::Float64 = 5.0,
	distance::Float64 = 0.4,
	power::Integer = -1)

	# Strategy
	if (strategy == N2Plain)
		jet_reconstruction = sequential_jet_reconstruct
        strategy_name = "N2Plain"
	elseif (strategy == N2TiledSoAGlobal)
		jet_reconstruction = tiled_jet_reconstruct_soa_global
        strategy_name = "N2TiledSoAGlobal"
	elseif (strategy == N2TiledSoATile)
		jet_reconstruction = tiled_jet_reconstruct_soa_tile
        strategy_name = "N2TiledSoATile"
    elseif (strategy == N2TiledLL)
		jet_reconstruction = tiled_jet_reconstruct_ll
        strategy_name = "N2TiledLL"
	else
		throw(ErrorException("Strategy not yet implemented"))
	end

    # Now run our jet reconstruction...
    events::Vector{Vector{PseudoJet}} = read_final_state_particles("test/data/events.hepmc3")
    if strategy == N2TiledLL
		event_vector = events
	else
		# First, convert all events into the Vector of Vectors that Atell's
		# code likes
		event_vector = pseudojets2vectors(events)
	end
    # event_vector = pseudojets2vectors(events)
    jet_collection = FinalJets[]
    for (ievt, event) in enumerate(event_vector)
        finaljets, _ = jet_reconstruction(event, R=distance, p=power)
        fj = final_jets(finaljets, ptmin)
        sort_jets!(fj)
        push!(jet_collection, FinalJets(ievt, fj))
    end

    @testset "Jet Reconstruction, $strategy_name" begin
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
                        @test jet.rap ≈ fastjet_jets[ievt]["jets"][ijet]["rap"] atol=1e-7
                        @test jet.phi ≈ fastjet_jets[ievt]["jets"][ijet]["phi"] atol=1e-7
                        @test jet.pt ≈ fastjet_jets[ievt]["jets"][ijet]["pt"] rtol=1e-6
                    end
                end
            end
        end
    end
end

"""Original test implementation"""
function original_tests()
    """
    A function for test comparison
    """
    function arrcompare(y, yt; eps=0.1)
        for i in eachindex(y)
            if !(sum(abs.(y[i] .- yt[i]) .< eps) == length(y[i]))
                return false
            end
        end
        true
    end

    NUMBER_OF_TESTS = 12 # number of test files in the data folder

    @testset "N2Plain Original Tests" begin
        # @test anti_kt(X) == Y

        # additional simple test
        @test arrcompare(
            anti_kt_algo([
                [99.0, 0.1, 0, 100],
                [4.0, -0.1, 0, 5.0],
                [-99., 0., 0.0, 99]
            ], R=0.7)[1],

            [[103.0, 0.0, 0.0, 105.0], [-99.0, 0.0, 0.0, 99.0]]
        )

        # actual testing
        for i in 1:NUMBER_OF_TESTS
            istr = string(i)

            @test arrcompare(
                    sort(
                        anti_kt_algo(
                            loadjets("test/data/"*istr*".dat"), R=1
                        )[1], lt=(a,b)->(pt(a)>pt(b))
                    ),

                    loadjets("test/data/"*istr*"-fj-result.dat")
                )
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
