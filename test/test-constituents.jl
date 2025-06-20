# Tests of jet constituent retrieval and parentage

include("common.jl")

# PseudoJet comparison test
function Base.isapprox(j1::PseudoJet, j2::PseudoJet)
    isapprox(j1.E, j2.E) && isapprox(j1.px, j2.px) &&
        isapprox(j1.py, j2.py) && isapprox(j1.pz, j2.pz)
end

# Expected constituent indexes and parent indexes
const expected_constituent_indexes = [84, 85, 139, 86, 133, 74, 79, 124, 76, 75, 163]
const expected_parent_indexes = [320, 335]

input_file = events_file_pp
events = read_final_state_particles(input_file)

# Event to pick
event_no = 1

cluster_seq = jet_reconstruct(events[event_no]; algorithm = JetAlgorithm.Kt, R = 1.0)

# Retrieve the exclusive pj_jets, but as `PseudoJet` types
pj_jets = inclusive_jets(cluster_seq; ptmin = 5.0, T = PseudoJet)

@testset "Jet constituents" begin
    @testset "Constituents of jet number $(event_no)" begin
        my_constituents = JetReconstruction.constituents(pj_jets[event_no], cluster_seq)
        @test size(my_constituents)[1] == 11
        for (i, idx) in enumerate(expected_constituent_indexes)
            @test my_constituents[i] ≈ events[1][idx]
        end
    end

    @testset "Constituent indexes for jet number $(event_no)" begin
        my_constituent_indexes = constituent_indexes(pj_jets[event_no], cluster_seq)
        @test size(my_constituent_indexes)[1] == 11
        # Testing the index values is sufficient, the content came from the original input file!
        @test my_constituent_indexes == expected_constituent_indexes
    end
end

@testset "Jet parents" begin
    @testset "Parent of jet number $(event_no)" begin
        my_parents = JetReconstruction.parent_jets(pj_jets[event_no], cluster_seq)
        @test my_parents[1] ≈
              cluster_seq.jets[cluster_seq.history[expected_parent_indexes[1]].jetp_index]
        @test my_parents[2] ≈
              cluster_seq.jets[cluster_seq.history[expected_parent_indexes[2]].jetp_index]
    end
    @testset "Parents of input cluster" begin
        no_parents = JetReconstruction.parent_jets(cluster_seq.jets[1], cluster_seq)
        @test isnothing(no_parents[1]) && isnothing(no_parents[2])
    end
end
