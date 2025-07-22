# Contains all common functions and necessary imports
include("common.jl")

events = JetReconstruction.read_final_state_particles(events_file_pp)

methods = Dict("jet-filter" => jet_filtering,
               "jet-trim" => jet_trimming,
               "mass-drop" => mass_drop,
               "soft-drop" => soft_drop)

for m in keys(methods)
    test_file = joinpath(@__DIR__, "data",
                         "fastjet-$m.json.zst")

    test_data = read_fastjet_outputs(test_file)
    tol = 1e-7

    if m === "jet-filter"
        groomer = (radius = 0.3, hardest_jets = 3)
        testname = "FastJet comparison: Jet Filtering, R=$(groomer.radius), n=$(groomer.hardest_jets)"

    elseif m === "jet-trim"
        groomer = (radius = 0.3, fraction = 0.3, recluster_method = JetAlgorithm.CA)
        testname = "FastJet comparison: Jet Trimming, R=$(groomer.radius), fraction=$(groomer.fraction), alg=$(groomer.recluster_method)"

    elseif m === "mass-drop"
        groomer = (mu = 0.67, y = 0.09)
        testname = "FastJet comparison: Mass Drop Tagging, μ=$(groomer.mu), y=$(groomer.y)"

    else
        groomer = (zcut = 0.1, beta = 2.0)
        testname = "FastJet comparison: Soft Drop Tagging, zcut=$(groomer.zcut), β=$(groomer.beta)"
        tol = 1e-4
    end

    @testset "$testname" begin
        for (ievt, evt) in enumerate(events)
            cluster_seq = jet_reconstruct(evt; algorithm = JetAlgorithm.CA, R = 1.0)
            jets = inclusive_jets(cluster_seq, PseudoJet; ptmin = 5.0)
            groomed = sort!([methods[m](jet, cluster_seq; groomer...)
                             for jet in jets], by = JetReconstruction.pt2, rev = true)

            @test length(groomed) === length(test_data[ievt]["jets"])
            for (ijet, jet) in enumerate(groomed)
                @test JetReconstruction.rapidity(jet)≈test_data[ievt]["jets"][ijet]["rap"] atol=tol
                @test JetReconstruction.phi(jet)≈test_data[ievt]["jets"][ijet]["phi"] atol=tol
                @test JetReconstruction.pt(jet)≈test_data[ievt]["jets"][ijet]["pt"] atol=tol
            end
        end
    end
end
