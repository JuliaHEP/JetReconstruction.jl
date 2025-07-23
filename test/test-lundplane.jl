# Contains all common functions and necessary imports
include("common.jl")

events = JetReconstruction.read_final_state_particles(events_file_pp)

test_file = joinpath(@__DIR__, "data",
                     "fastjet-lundplane.json.zst")

test_data = read_fastjet_outputs(test_file)
tol = 1e-4
count = 1

@testset "FastJet comparison: Lund Plane coordinates" begin
    for (ievt, evt) in enumerate(events)
        cluster_seq = jet_reconstruct(evt; algorithm = JetAlgorithm.AntiKt, R = 1.0)
        jets = sort!(inclusive_jets(cluster_seq; ptmin = 10.0, T = PseudoJet),
                     by = JetReconstruction.pt2, rev = true)

        for (ijet, jet) in enumerate(jets)
            lundvars = generate_lund_emissions(jet, cluster_seq)
            for temp in lundvars
                @test temp.kt≈test_data[count]["kt"] atol=tol
                @test temp.delta≈test_data[count]["Delta"] atol=tol
                @test temp.h_pt≈test_data[count]["h_pt"] atol=tol
                @test temp.s_pt≈test_data[count]["s_pt"] atol=tol
                global count += 1
            end
        end
    end
end
