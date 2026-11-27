include("common.jl")

using Random

@testset "FastJet explicit ghost reference" begin
    reference = JSON.parse(read(JetReconstruction.open_with_stream(joinpath(@__DIR__,
                                                                            "data",
                                                                            "jet-areas-fastjet.json.zst")),
                                String))
    input = vcat(reference["particles"], reference["ghosts"])
    particles = [PseudoJet(momentum...; cluster_hist_index = index)
                 for (index, momentum) in enumerate(input)]
    for result in reference["results"],
        strategy in (RecoStrategy.N2Plain, RecoStrategy.N2Tiled)

        algorithm = getproperty(JetAlgorithm, Symbol(result["algorithm"]))
        recombine = (addjets_escheme, addjets_ptscheme, addjets_pt2scheme)[result["scheme"] + 1]
        cs = jet_reconstruct(particles; algorithm, p = 1.5, R = reference["R"],
                             strategy, recombine, preprocess = nothing)
        data = JetReconstruction._postprocess_area(cs, length(reference["particles"]),
                                                   reference["ghost_area"]; recombine)
        csa = ClusterSequenceArea(cs, data)
        jets = sort(inclusive_jets(csa; ptmin = 1); by = JetReconstruction.pt, rev = true)
        @test length(jets) == length(result["jets"])
        for (jet, expected) in zip(jets, result["jets"])
            @test [
                JetReconstruction.px(jet),
                JetReconstruction.py(jet),
                JetReconstruction.pz(jet),
                JetReconstruction.energy(jet)
            ] ≈ expected["momentum"]
            @test area(jet, csa) ≈ expected["area"]
            vector = area_4vector(jet, csa)
            @test [
                JetReconstruction.px(vector),
                JetReconstruction.py(vector),
                JetReconstruction.pz(vector),
                JetReconstruction.energy(vector)
            ] ≈
                  expected["area_4vector"]
        end
    end
end

@testset "active jet areas" begin
    particles = [PseudoJet(; pt = 100.0, rap = 0.0, phi = 0.0),
        PseudoJet(; pt = 80.0, rap = 0.1, phi = π)]
    original_particles = copy(particles)
    spec = GhostedAreaSpec(ghost_maxrap = 1.0, ghost_area = 0.02, seed = 12345)

    @testset "algorithms and strategies" begin
        for algorithm in (JetAlgorithm.AntiKt, JetAlgorithm.CA, JetAlgorithm.Kt,
                          JetAlgorithm.GenKt)
            results = map((RecoStrategy.N2Plain, RecoStrategy.N2Tiled, RecoStrategy.Best)) do strategy
                csa = jet_reconstruct_area(particles; algorithm, p = 1.5, R = 0.4,
                                           area_spec = spec, strategy)
                jets = sort(inclusive_jets(csa; ptmin = 1.0); by = JetReconstruction.pt,
                            rev = true)
                @test jets isa Vector{PseudoJet}
                @test length(jets) == 2
                @test all(!is_pure_ghost(jet, csa) for jet in jets)
                @test all(area(jet, csa) > 0 for jet in jets)
                @test all(area_error(jet, csa) == 0 for jet in jets)
                @test total_area(csa) ≈ 4π * spec.ghost_maxrap
                @test n_hard_particles(csa) == length(particles)
                @test has_explicit_ghosts(csa)
                @test !has_dangerous_particles(csa)
                @test max_ghost_pt2(csa) > 0
                all_jets = inclusive_jets(csa; include_pure_ghosts = true)
                @test length(all_jets) > length(jets)
                @test sum(area(jet, csa) for jet in all_jets) ≈ total_area(csa)
                if algorithm == JetAlgorithm.AntiKt
                    @test all(isapprox(area(jet, csa), π * 0.4^2; rtol = 0.2)
                              for jet in jets)
                end
                [area(jet, csa) for jet in jets]
            end
            @test results[1] ≈ results[2]
            @test results[1] ≈ results[3]
        end
        @test particles == original_particles
    end

    @testset "grid and randomness" begin
        grid = JetReconstruction.GhostGridSpec(spec)
        @test grid.n_ghosts == 2 * grid.nrap * grid.nphi
        @test grid.n_ghosts * grid.actual_ghost_area ≈ 4π * spec.ghost_maxrap
        ghosts, ghost_area = JetReconstruction._make_area_ghosts(spec, MersenneTwister(7))
        @test length(ghosts) == grid.n_ghosts
        @test ghost_area == grid.actual_ghost_area
        @test all(-spec.ghost_maxrap <= JetReconstruction.rapidity(ghost) <=
                  spec.ghost_maxrap
                  for ghost in ghosts)
        @test all(0 <= JetReconstruction.phi(ghost) < 2π for ghost in ghosts)
        @test all(isapprox(JetReconstruction.energy(ghost),
                           JetReconstruction.pt(ghost) *
                           cosh(JetReconstruction.rapidity(ghost)))
                  for ghost in ghosts)
        regular = GhostedAreaSpec(ghost_maxrap = 1.0, ghost_area = 0.02,
                                  grid_scatter = 0, kt_scatter = 0)
        regular_ghosts, _ = JetReconstruction._make_area_ghosts(regular, MersenneTwister(1))
        @test all(JetReconstruction.pt(ghost) ≈ regular.mean_ghost_pt
                  for ghost in regular_ghosts)
        @test JetReconstruction.rapidity(first(regular_ghosts)) ≈ -1 + grid.drap / 2
        @test JetReconstruction.phi(first(regular_ghosts)) ≈ grid.dphi / 2
        first_run = jet_reconstruct_area(particles; algorithm = JetAlgorithm.AntiKt,
                                         area_spec = spec)
        second_run = jet_reconstruct_area(particles; algorithm = JetAlgorithm.AntiKt,
                                          area_spec = spec)
        @test cluster_sequence(first_run).jets == cluster_sequence(second_run).jets
        override_run = jet_reconstruct_area(particles; algorithm = JetAlgorithm.AntiKt,
                                            area_spec = spec, rng = MersenneTwister(8))
        @test cluster_sequence(first_run).jets != cluster_sequence(override_run).jets
        unseeded = GhostedAreaSpec(ghost_maxrap = 1, ghost_area = 0.1)
        explicit_rng = jet_reconstruct_area(particles; algorithm = JetAlgorithm.CA,
                                            area_spec = unseeded, rng = MersenneTwister(8))
        repeated_rng = jet_reconstruct_area(particles; algorithm = JetAlgorithm.CA,
                                            area_spec = unseeded, rng = MersenneTwister(8))
        @test cluster_sequence(explicit_rng).jets == cluster_sequence(repeated_rng).jets
    end

    @testset "area vectors and ghost classification" begin
        csa = jet_reconstruct_area(particles; algorithm = JetAlgorithm.AntiKt,
                                   R = 0.4, area_spec = spec)
        cs = cluster_sequence(csa)
        data = csa.area_data
        for jet in inclusive_jets(csa; include_pure_ghosts = true)
            indexes = constituent_indexes(jet, cs)
            ghost_indexes = filter(index -> index > length(particles), indexes)
            @test is_pure_ghost(jet, csa) == (length(ghost_indexes) == length(indexes))
            @test area(jet, csa) ≈ length(ghost_indexes) * data.ghost_area
            components = map((JetReconstruction.px, JetReconstruction.py,
                              JetReconstruction.pz, JetReconstruction.energy)) do component
                sum(component(cs.jets[index]) * data.ghost_area /
                    JetReconstruction.pt(cs.jets[index])
                    for index in ghost_indexes; init = 0.0)
            end
            vector = area_4vector(jet, csa)
            @test [
                JetReconstruction.px(vector),
                JetReconstruction.py(vector),
                JetReconstruction.pz(vector),
                JetReconstruction.energy(vector)
            ] ≈ collect(components)
        end
        ghost_only = jet_reconstruct_area(PseudoJet[]; algorithm = JetAlgorithm.Kt,
                                          area_spec = spec)
        @test isempty(inclusive_jets(ghost_only))
        @test n_hard_particles(ghost_only) == 0
        @test all(is_pure_ghost(jet, ghost_only)
                  for jet in inclusive_jets(ghost_only; include_pure_ghosts = true))
        dangerous = jet_reconstruct_area(particles; algorithm = JetAlgorithm.AntiKt,
                                         area_spec = GhostedAreaSpec(ghost_maxrap = 1,
                                                                     ghost_area = 0.1,
                                                                     mean_ghost_pt = 1e-5))
        @test has_dangerous_particles(dangerous)
        for recombine in (addjets_escheme, addjets_ptscheme, addjets_pt2scheme)
            outside = [PseudoJet(; pt = 100, rap = 3, phi = 0),
                PseudoJet(; pt = 10, rap = 3.05, phi = 0)]
            outside_csa = jet_reconstruct_area(outside; algorithm = JetAlgorithm.AntiKt,
                                               R = 0.4, area_spec = spec, recombine)
            jet = only(inclusive_jets(outside_csa; ptmin = 1))
            @test area(jet, outside_csa) == 0
            vector = area_4vector(jet, outside_csa)
            @test [
                JetReconstruction.px(vector),
                JetReconstruction.py(vector),
                JetReconstruction.pz(vector),
                JetReconstruction.energy(vector)
            ] == zeros(4)
            inside_csa = jet_reconstruct_area(particles; algorithm = JetAlgorithm.AntiKt,
                                              R = 0.4, area_spec = spec, recombine)
            @test all(isfinite(JetReconstruction.energy(area_4vector(jet, inside_csa)))
                      for jet in inclusive_jets(inside_csa))
        end
    end

    @testset "inputs and output conversions" begin
        raw = jet_reconstruct_area(particles; algorithm = JetAlgorithm.AntiKt,
                                   area_spec = spec, preprocess = nothing)
        converted = jet_reconstruct_area(lorentzvector.(particles);
                                         algorithm = JetAlgorithm.AntiKt, area_spec = spec)
        raw_jets = inclusive_jets(raw; ptmin = 1)
        converted_jets = inclusive_jets(converted; ptmin = 1)
        @test all(component.(raw_jets) ≈ component.(converted_jets)
                  for component in (JetReconstruction.px, JetReconstruction.py,
                                    JetReconstruction.pz, JetReconstruction.energy))
        @test area.(raw_jets, Ref(raw)) ≈ area.(converted_jets, Ref(converted))
        for output_type in (LorentzVector, LorentzVectorCyl,
                            LorentzVector{Float32}, LorentzVectorCyl{Float32})
            jets = inclusive_jets(raw, output_type; ptmin = 1)
            @test eltype(jets) == output_type
            @test length(jets) == 2
            @test all(JetReconstruction.pt(jet) > 1 for jet in jets)
            @test eltype(inclusive_jets(raw, output_type; ptmin = 1e10)) == output_type
        end
        @test_throws ArgumentError inclusive_jets(raw, Int)
        foreign = PseudoJet(; pt = 1, rap = 0, phi = 0, cluster_hist_index = 1)
        for query in (area, area_4vector, area_error, is_pure_ghost)
            @test_throws ArgumentError query(first(particles), raw)
            @test_throws ArgumentError query(foreign, raw)
        end
    end

    @testset "invalid configurations" begin
        for name in (:ghost_area, :ghost_maxrap, :mean_ghost_pt), value in (0, -1, Inf, NaN)
            invalid_spec = GhostedAreaSpec(; NamedTuple{(name,)}((value,))...)
            @test_throws ArgumentError JetReconstruction.GhostGridSpec(invalid_spec)
        end
        for name in (:grid_scatter, :kt_scatter), value in (-1, Inf, NaN)
            invalid_spec = GhostedAreaSpec(; NamedTuple{(name,)}((value,))...)
            @test_throws ArgumentError JetReconstruction.GhostGridSpec(invalid_spec)
        end
        for invalid_spec in (GhostedAreaSpec(repeat = 0), GhostedAreaSpec(repeat = 2),
                             GhostedAreaSpec(kt_scatter = 2),
                             GhostedAreaSpec(mean_ghost_pt = 1e-200),
                             GhostedAreaSpec(mean_ghost_pt = 1e200),
                             GhostedAreaSpec(ghost_area = 1e-100),
                             GhostedAreaSpec(ghost_area = 1e-20),
                             GhostedAreaSpec(ghost_maxrap = 1000))
            @test_throws ArgumentError JetReconstruction.GhostGridSpec(invalid_spec)
        end
        for radius in (0, -1, Inf, NaN, 1e200, 1e-200)
            @test_throws ArgumentError jet_reconstruct_area(particles;
                                                            algorithm = JetAlgorithm.AntiKt,
                                                            R = radius, area_spec = spec)
        end
        for algorithm in (JetAlgorithm.EEKt, JetAlgorithm.Durham, JetAlgorithm.Valencia)
            @test_throws ArgumentError jet_reconstruct_area(particles; algorithm,
                                                            area_spec = spec)
        end
        for power in (-2.0, 2.0)
            @test_throws ArgumentError jet_reconstruct_area(particles;
                                                            algorithm = JetAlgorithm.GenKt,
                                                            p = power, area_spec = spec)
        end
        @test_throws ArgumentError jet_reconstruct_area(particles;
                                                        algorithm = JetAlgorithm.GenKt,
                                                        area_spec = spec)
        @test_throws ArgumentError jet_reconstruct_area(particles;
                                                        algorithm = JetAlgorithm.GenKt,
                                                        p = Inf, area_spec = spec)
    end

    @testset "azimuthal tile boundary" begin
        tiling = JetReconstruction.setup_tiling([-1.0, 1.0], 0.4)
        for angle in (0.0, prevfloat(Float64(2π)), Float64(2π), nextfloat(Float64(2π))),
            rap in (-10.0, 0.0, 10.0)

            index = JetReconstruction.tile_index(tiling, rap, angle)
            @test 1 <= index <= tiling._n_tiles
        end
        boundary_particle = PseudoJet(100.0, -1e-15, 0.0, 100.0)
        @test JetReconstruction.phi(boundary_particle) == 2π
        boundary_csa = jet_reconstruct_area([boundary_particle];
                                            algorithm = JetAlgorithm.AntiKt, R = 0.4,
                                            area_spec = spec,
                                            strategy = RecoStrategy.N2Tiled)
        @test length(inclusive_jets(boundary_csa; ptmin = 1)) == 1
    end
end
