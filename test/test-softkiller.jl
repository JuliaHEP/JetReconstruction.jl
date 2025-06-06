# Tests for SoftKiller algorithm

include("common.jl")

using Random
import JetReconstruction: SoftKiller, tile_index, softkiller_apply, select_ABS_RAP_max, PseudoJet

@testset "SoftKiller grid setup" begin
    # Test default constructor
    sk = SoftKiller(-2.5, 2.5, 0.4, 0.4)
    @test sk._ymin == -2.5
    @test sk._ymax == 2.5
    @test sk._requested_drap == 0.4
    @test sk._requested_dphi == 0.4
    @test sk._ntotal > 0

    # Test constructor with only ymin and square grid size
    sk2 = SoftKiller(2.5, 0.4)
    @test sk2._ymin == -2.5
    @test sk2._ymax == 2.5
    @test sk2._requested_drap == 0.4
    @test sk2._requested_dphi == 0.4
end

@testset "SoftKiller tile_index and ABS_RAP_max" begin
    sk = SoftKiller(-2.5, 2.5, 0.4, 0.4)

    # Test tile_index with a pseudojet within bounds
    pj1 = PseudoJet(pt=1.0, rap=0.0, phi=0.0, m=0.0)
    idx = tile_index(sk, pj1)
    @test (1 <= idx) && (idx <= sk._ntotal)
    
    # Test tile_index with a pseudojet beyond bounds
    pj2 = PseudoJet(pt=1.0, rap=10.0, phi=0.0, m=0.0)
    idx2 = tile_index(sk, pj2)
    @test idx2 == -1

    event = [pj1, pj2]

    # Test with an eta_max that should filter out pj2
    filtered = select_ABS_RAP_max(event, 2.0)
    @test pj1 in filtered
    @test !(pj2 in filtered)

    # Test with an eta_max that should keep both
    filtered2 = select_ABS_RAP_max(event, 10.0)
    @test pj1 in filtered2
    @test pj2 in filtered2
end

@testset "SoftKiller apply" begin
    N = 200
    rapmin, rapmax = -2.5, 2.5
    rng = Random.default_rng()
    Random.seed!(rng, 12535)
    sk = SoftKiller(rapmin, rapmax, 0.4, 0.4)
    # Exponentially distributed pt, uniform in rapidity and phi
    pts = randexp(rng, N) .+ 0.1
    raps = rand(rng, N) .* (rapmax - rapmin) .+ rapmin
    phis = rand(rng, N) .* (2π)

    # Helper function to create a PseudoJet from pt, eta, and phi
    function make_pj(pt, rap, phi)
        return PseudoJet(pt=pt, rap=rap, phi=phi, m=0.0)
    end
    event = [make_pj(pts[i], raps[i], phis[i]) for i in 1:N]

    reduced, threshold = softkiller!(sk, event)

    @test threshold > 0
    @test threshold < maximum(pts)
    @test length(reduced) < N

    # Check that only jets with pt > threshold remain
    for (pj, pt) in zip(event, pts)
        if pt > threshold
            @test pj in reduced
        else
            @test pj ∉ reduced
        end
    end
end

