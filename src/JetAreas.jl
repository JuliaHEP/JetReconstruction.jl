"""
    GhostedAreaSpec(; kwargs...)

Configuration for active jet-area calculation using explicit ghost particles.

The ghost grid covers rapidities from approximately `-ghost_maxrap` to
`+ghost_maxrap`. The scalar area represented by each ghost is controlled by
`ghost_area`.

This first implementation supports the FastJet-style active area with explicit
ghosts and `repeat == 1`.
"""

#These defaults are consistent with commonly used FastJet active-area settings
Base.@kwdef struct GhostedAreaSpec
    ghost_maxrap::Float64 = 6.0
    repeat::Int = 1
    ghost_area::Float64 = 0.01
    grid_scatter::Float64 = 1.0
    kt_scatter::Float64 = 0.1
    mean_ghost_pt::Float64 = 1.0e-100
    seed::Union{Nothing,Int} = nothing

end

struct GhostGridSpec
    nphi::Int
    nrap::Int
    dphi::Float64
    drap::Float64
    actual_ghost_area::Float64
    n_ghosts::Int
end

function _ghost_grid_spec(spec::GhostedAreaSpec) :: GhostGridSpec
    spec.repeat == 1 || throw(ArgumentError("jet area repeat > 1 is not implemented yet"))

    requested_spacing = sqrt(spec.ghost_area)

    nphi = max(1, floor(Int, 2π / requested_spacing + 0.5))
    dphi = 2π / nphi

    nrap = max(1, floor(Int, spec.ghost_maxrap / requested_spacing + 0.5))
    drap = spec.ghost_maxrap / nrap
    #FastJet does not use the requested ghost area exactly. 
    #It uses it to choose a grid spacing, then recomputes the actual ghost area from the final rapidity and phi spacing.
    actual_ghost_area = dphi * drap
    n_ghosts = 2 * nrap * nphi

    return GhostGridSpec(nphi, nrap, dphi, drap, actual_ghost_area, n_ghosts)
end

_copy_as_pseudojet(jet; cluster_hist_index::Int) = PseudoJet(
    px(jet),
    py(jet),
    pz(jet),
    E(jet);
    cluster_hist_index = cluster_hist_index,
)

_zero_pseudojet(; cluster_hist_index::Int = 0) = PseudoJet(
    0.0,
    0.0,
    0.0,
    0.0;
    cluster_hist_index = cluster_hist_index,
)

function _scale_as_pseudojet(jet, scale::Real; cluster_hist_index::Int)
    return PseudoJet(
        scale * px(jet),
        scale * py(jet),
        scale * pz(jet),
        scale * E(jet);
        cluster_hist_index = cluster_hist_index,
    )
end

function _make_area_ghosts(spec::GhostedAreaSpec)
    grid = _ghost_grid_spec(spec)
    rng = isnothing(spec.seed) ? Random.default_rng() : Random.MersenneTwister(spec.seed)

    ghosts = Vector{PseudoJet}()
    sizehint!(ghosts, grid.n_ghosts)

    #FastJet-style FJ3 layout: irap runs from -nrap to nrap - 1.
    for irap in -grid.nrap:(grid.nrap - 1)
        for iphi in 0:(grid.nphi - 1)
            phi = (iphi + 0.5) * grid.dphi
            phi += grid.dphi * (rand(rng) - 0.5) * spec.grid_scatter
            phi = mod(phi, 2π)

            rap = (irap + 0.5) * grid.drap
            rap += grid.drap * (rand(rng) - 0.5) * spec.grid_scatter

            kt = spec.mean_ghost_pt * (1 + (rand(rng) - 0.5) * spec.kt_scatter)

            ghost_px = kt * cos(phi)
            ghost_py = kt * sin(phi)

            pplus = kt * exp(rap)
            pminus = kt * exp(-rap)

            ghost_pz = 0.5 * (pplus - pminus)
            ghost_e = 0.5 * (pplus + pminus)

            push!(
                ghosts,
                PseudoJet(
                    ghost_px,
                    ghost_py,
                    ghost_pz,
                    ghost_e;
                    cluster_hist_index = 0,
                ),
            )
        end
    end

    return ghosts, grid.actual_ghost_area
end

struct JetAreaData
    ghost_area::Float64
    n_ghosts::Int
    n_hard_particles::Int
    areas::Vector{Float64}
    area_4vectors::Vector{PseudoJet}
    is_pure_ghost::Vector{Bool}
    max_ghost_pt2::Float64
    has_dangerous_particles::Bool
end

struct ClusterSequenceArea
    clusterseq::ClusterSequence{PseudoJet}
    area_data::JetAreaData
end

function _postprocess_area(
    cs::ClusterSequence{PseudoJet},
    n_hard_particles::Int,
    ghost_area::Float64;
    recombine = addjets_escheme,
)
    n_initial = cs.n_initial_jets
    n_history = length(cs.history)

    areas = zeros(Float64, n_history)
    area_4vectors = Vector{PseudoJet}(undef, n_history)
    is_pure_ghost = falses(n_history)

    for i in 1:n_history
        area_4vectors[i] = _zero_pseudojet(cluster_hist_index = i)
    end

    n_ghosts = n_initial - n_hard_particles

    max_ghost_pt2 = 0.0

    for i in 1:n_initial
        jet = cs.jets[i]

        if i > n_hard_particles
            is_pure_ghost[i] = true
            areas[i] = ghost_area

            jet_pt = pt(jet)
            area_4vectors[i] = _scale_as_pseudojet(
                jet,
                ghost_area / jet_pt;
                cluster_hist_index = i,
            )

            max_ghost_pt2 = max(max_ghost_pt2, pt2(jet))
        else
            is_pure_ghost[i] = false
            areas[i] = 0.0
            area_4vectors[i] = _zero_pseudojet(cluster_hist_index = i)
        end
    end

    has_dangerous_particles = false
    if max_ghost_pt2 > 0.0
        danger_ratio = eps(Float64)^2

        for i in 1:n_hard_particles
            if danger_ratio * pt2(cs.jets[i]) <= max_ghost_pt2
                has_dangerous_particles = true
                break
            end
        end
    end

    for i in (n_initial + 1):n_history
        history_element = cs.history[i]
        parent1 = history_element.parent1
        parent2 = history_element.parent2

        if parent2 == BeamJet
            is_pure_ghost[i] = is_pure_ghost[parent1]
            areas[i] = areas[parent1]
            area_4vectors[i] = area_4vectors[parent1]
        elseif parent2 > 0
            is_pure_ghost[i] = is_pure_ghost[parent1] && is_pure_ghost[parent2]
            areas[i] = areas[parent1] + areas[parent2]
            area_4vectors[i] = recombine(
                area_4vectors[parent1],
                area_4vectors[parent2];
                cluster_hist_index = i,
            )
        else
            areas[i] = 0.0
            area_4vectors[i] = _zero_pseudojet(cluster_hist_index = i)
            is_pure_ghost[i] = false
        end
    end

    return JetAreaData(
        ghost_area,
        n_ghosts,
        n_hard_particles,
        areas,
        area_4vectors,
        is_pure_ghost,
        max_ghost_pt2,
        has_dangerous_particles,
    )
end

"""
    jet_reconstruct_area(particles; algorithm, R, area_spec = GhostedAreaSpec(), kwargs...)

Cluster particles with explicit ghosts and compute active jet areas.

This currently supports pp-style algorithms only.
"""
function jet_reconstruct_area(
    particles;
    algorithm::JetAlgorithm.Algorithm,
    R = 1.0,
    p = nothing,
    area_spec::GhostedAreaSpec = GhostedAreaSpec(),
    strategy = RecoStrategy.N2Plain,
    recombine = addjets_escheme,
    preprocess = preprocess_escheme,
)
    if algorithm ∉ (
        JetAlgorithm.Kt,
        JetAlgorithm.CA,
        JetAlgorithm.AntiKt,
        JetAlgorithm.GenKt,
    )
        throw(ArgumentError("jet_reconstruct_area currently supports pp algorithms only"))
    end
    if strategy != RecoStrategy.N2Plain
        throw(ArgumentError("jet_reconstruct_area currently supports only RecoStrategy.N2Plain"))
    end

    hard_particles = Vector{PseudoJet}()
    sizehint!(hard_particles, length(particles))

    for (i, particle) in enumerate(particles)
        hard_particle =
            isnothing(preprocess) ?
            _copy_as_pseudojet(particle; cluster_hist_index = i) :
            preprocess(particle, PseudoJet; cluster_hist_index = i)

        push!(hard_particles, hard_particle)
    end

    ghosts, actual_ghost_area = _make_area_ghosts(area_spec)

    n_hard_particles = length(hard_particles)

    for i in eachindex(ghosts)
        ghosts[i] = _copy_as_pseudojet(
            ghosts[i];
            cluster_hist_index = n_hard_particles + i,
        )
    end

    particles_with_ghosts = vcat(hard_particles, ghosts)

    cs = jet_reconstruct(
        particles_with_ghosts;
        algorithm = algorithm,
        R = R,
        p = p,
        strategy = strategy,
        recombine = recombine,
        preprocess = nothing,
    )

    area_data = _postprocess_area(
        cs,
        n_hard_particles,
        actual_ghost_area;
        recombine = recombine,
    )

    return ClusterSequenceArea(cs, area_data)
end

cluster_sequence(csa::ClusterSequenceArea) = csa.clusterseq

has_explicit_ghosts(::ClusterSequenceArea) = true

function area(jet, csa::ClusterSequenceArea)
    return csa.area_data.areas[cluster_hist_index(jet)]
end

function area_4vector(jet, csa::ClusterSequenceArea)
    return csa.area_data.area_4vectors[cluster_hist_index(jet)]
end

function area_error(jet, csa::ClusterSequenceArea)
    return 0.0
end

function is_pure_ghost(jet, csa::ClusterSequenceArea)
    return csa.area_data.is_pure_ghost[cluster_hist_index(jet)]
end

function total_area(csa::ClusterSequenceArea)
    return csa.area_data.n_ghosts * csa.area_data.ghost_area
end

function max_ghost_pt2(csa::ClusterSequenceArea)
    return csa.area_data.max_ghost_pt2
end

function has_dangerous_particles(csa::ClusterSequenceArea)
    return csa.area_data.has_dangerous_particles
end

function n_hard_particles(csa::ClusterSequenceArea)
    return csa.area_data.n_hard_particles
end

function _convert_area_result_jet(jet::PseudoJet, ::Type{PseudoJet})
    return jet
end

function _convert_area_result_jet(jet::PseudoJet, ::Type{T}) where {T <: LorentzVector}
    return lorentzvector(jet)
end

function _convert_area_result_jet(jet::PseudoJet, ::Type{T}) where {T <: LorentzVectorCyl}
    return lorentzvector_cyl(jet)
end

function inclusive_jets(
    csa::ClusterSequenceArea,
    ::Type{T} = PseudoJet;
    ptmin = 0.0,
    include_pure_ghosts::Bool = false,
) where {T}
    jets = inclusive_jets(csa.clusterseq, PseudoJet; ptmin = ptmin)

    if !include_pure_ghosts
        jets = filter(jet -> !is_pure_ghost(jet, csa), jets)
    end

    return [_convert_area_result_jet(jet, T) for jet in jets]
end
