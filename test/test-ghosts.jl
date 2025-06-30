include("common.jl")

@testset "Ghost Creation Tests" begin
    resolution = 1000
    gen = GhostedArea(resolution)

    # Create original event with preexisting PseudoJets
    event = [
        PseudoJet(1.1, 1.0, 1.3, 1.2),
        PseudoJet(2.4, 3.2, 0.8, 1.7),
        PseudoJet(0.9, 0.5, 3.0, 4.3),
        PseudoJet(0.9, 6.0, 5.2, 0.4)
    ]
    num_original = length(event)

    # Use the function being tested
    add_ghosts!(event, gen)

    # Ensure that the correct number of ghosts were generated
    num_ghosts = resolution * resolution
    @test length(event) == num_original + num_ghosts

    # Ensure that the original jets have not changed
    @test event[1].px == 1.1
    @test event[2].py == 3.2
    @test event[3].pz == 3.0
    @test event[4].E == 0.4

    # Check that all added ghosts have the correct pt, use an error bound due to floating point error
    @test all(i -> i._pt2 - 1.0e-90 <= 1.0e-89, event[(num_original + 1):end])
end
