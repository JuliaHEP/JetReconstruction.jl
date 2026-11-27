"""
    GhostedAreaSpec(; kwargs...)

Configuration for active jet-area calculation using explicit ghost particles.

The ghost grid covers rapidities from approximately `-ghost_maxrap` to
`+ghost_maxrap`. Choose coverage at least `R` beyond the jets of interest to
avoid acceptance edge effects. The requested `ghost_area` determines the grid
spacing; the actual area is adjusted to fit the rapidity and azimuthal ranges.

Keyword defaults are `ghost_maxrap = 6.0`, `ghost_area = 0.01`,
`grid_scatter = 1.0` (position scatter relative to a grid cell),
`kt_scatter = 0.1` (relative transverse momentum scatter), and
`mean_ghost_pt = 1e-100`. All must be finite, the coverage and area positive,
`grid_scatter ≥ 0`, and `0 ≤ kt_scatter < 2`. Ghost transverse momenta must
have finite squares greater than `1e-300` for the clustering distance metric.

Only `repeat = 1` is supported for explicit ghosts. An optional integer `seed`
creates a local `MersenneTwister`; alternatively pass `rng` to
[`jet_reconstruct_area`](@ref). An explicit `rng` takes precedence over `seed`.
"""
Base.@kwdef struct GhostedAreaSpec
    ghost_maxrap::Float64 = 6.0
    repeat::Int = 1
    ghost_area::Float64 = 0.01
    grid_scatter::Float64 = 1.0
    kt_scatter::Float64 = 0.1
    mean_ghost_pt::Float64 = 1.0e-100
    seed::Union{Nothing, Int} = nothing
end

"""Grid dimensions and actual cell area for an explicit ghost distribution."""
struct GhostGridSpec
    nphi::Int
    nrap::Int
    dphi::Float64
    drap::Float64
    actual_ghost_area::Float64
    n_ghosts::Int
end

function GhostGridSpec(spec::GhostedAreaSpec)
    spec.repeat == 1 || throw(ArgumentError("explicit ghost areas require repeat == 1"))
    for (name, value) in ((:ghost_maxrap, spec.ghost_maxrap),
                          (:ghost_area, spec.ghost_area),
                          (:mean_ghost_pt, spec.mean_ghost_pt))
        isfinite(value) && value > 0 ||
            throw(ArgumentError("$name must be finite and positive"))
    end
    isfinite(spec.grid_scatter) && spec.grid_scatter >= 0 ||
        throw(ArgumentError("grid_scatter must be finite and nonnegative"))
    isfinite(spec.kt_scatter) && 0 <= spec.kt_scatter < 2 ||
        throw(ArgumentError("kt_scatter must be finite and in [0, 2)"))
    min_pt2 = (spec.mean_ghost_pt * (1 - spec.kt_scatter / 2))^2
    max_pt2 = (spec.mean_ghost_pt * (1 + spec.kt_scatter / 2))^2
    min_pt2 > 1e-300 && isfinite(max_pt2) ||
        throw(ArgumentError("ghost transverse momentum squares must be finite and greater than 1e-300"))

    requested_spacing = sqrt(spec.ghost_area)

    phi_cells = max(1.0, floor(2π / requested_spacing + 0.5))
    rap_cells = max(1.0, floor(spec.ghost_maxrap / requested_spacing + 0.5))
    phi_cells < typemax(Int) && rap_cells < typemax(Int) ||
        throw(ArgumentError("ghost grid dimensions exceed the integer range"))
    nphi = Int(phi_cells)
    dphi = 2π / nphi

    nrap = Int(rap_cells)
    drap = spec.ghost_maxrap / nrap
    actual_ghost_area = dphi * drap
    nrap <= div(typemax(Int), 2) ÷ nphi ||
        throw(ArgumentError("ghost count exceeds the integer range"))
    n_ghosts = 2 * nrap * nphi
    max_rap = spec.ghost_maxrap + max(0.0, spec.grid_scatter - 1) * drap / 2
    max_energy = spec.mean_ghost_pt * (1 + spec.kt_scatter / 2) * cosh(max_rap)
    isfinite(max_energy^2) && isfinite(actual_ghost_area) && actual_ghost_area > 0 ||
        throw(ArgumentError("ghost grid produces nonfinite momenta or cell areas"))

    return GhostGridSpec(nphi, nrap, dphi, drap, actual_ghost_area, n_ghosts)
end

function _zero_pseudojet(; cluster_hist_index::Int = 0)
    PseudoJet(0.0,
              0.0,
              0.0,
              0.0;
              cluster_hist_index = cluster_hist_index,)
end

function _scale_as_pseudojet(jet, scale::Real; cluster_hist_index::Int)
    return PseudoJet(scale * px(jet),
                     scale * py(jet),
                     scale * pz(jet),
                     scale * E(jet);
                     cluster_hist_index = cluster_hist_index,)
end

function _rng_from_area_spec(spec::GhostedAreaSpec)
    isnothing(spec.seed) ? Random.default_rng() : Random.MersenneTwister(spec.seed)
end

"""Generate massless explicit ghosts, returning their actual scalar cell area."""
function _make_area_ghosts(spec::GhostedAreaSpec, rng::Random.AbstractRNG)
    grid = GhostGridSpec(spec)
    ghosts = Vector{PseudoJet}()
    sizehint!(ghosts, grid.n_ghosts)

    for irap in (-grid.nrap):(grid.nrap - 1)
        for iphi in 0:(grid.nphi - 1)
            phi = (iphi + 0.5) * grid.dphi
            phi += grid.dphi * (rand(rng) - 0.5) * spec.grid_scatter
            phi = mod(phi, 2π)

            rap = (irap + 0.5) * grid.drap
            rap += grid.drap * (rand(rng) - 0.5) * spec.grid_scatter

            kt = spec.mean_ghost_pt * (1 + (rand(rng) - 0.5) * spec.kt_scatter)

            push!(ghosts, PseudoJet(; pt = kt, rap = rap, phi = phi))
        end
    end

    return ghosts, grid.actual_ghost_area
end

"""
    JetAreaData

Explicit ghost metadata and scalar areas, area four-vectors, and pure-ghost flags
indexed by clustering history. Usually obtained through [`jet_reconstruct_area`](@ref).
"""
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

"""
    ClusterSequenceArea

A clustering sequence with explicit ghost area data. Use [`inclusive_jets`](@ref)
to obtain area-queryable `PseudoJet`s and [`cluster_sequence`](@ref) to access
the underlying history and constituents (which include ghosts).
"""
struct ClusterSequenceArea
    clusterseq::ClusterSequence{PseudoJet}
    area_data::JetAreaData
end

"""Propagate initial ghost areas and recombined area four-vectors through a history."""
function _postprocess_area(cs::ClusterSequence{PseudoJet},
                           n_hard_particles::Int,
                           ghost_area::Float64;
                           recombine = addjets_escheme,)
    n_initial = cs.n_initial_jets
    n_history = length(cs.history)

    areas = zeros(Float64, n_history)
    area_4vectors = fill(_zero_pseudojet(), n_history)
    is_pure_ghost = fill(false, n_history)

    n_ghosts = n_initial - n_hard_particles

    max_ghost_pt2 = 0.0

    for i in 1:n_initial
        jet = cs.jets[i]

        if i > n_hard_particles
            is_pure_ghost[i] = true
            areas[i] = ghost_area

            jet_pt = pt(jet)
            area_4vectors[i] = _scale_as_pseudojet(jet,
                                                   ghost_area / jet_pt;
                                                   cluster_hist_index = i,)

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
            if iszero(areas[parent1])
                area_4vectors[i] = PseudoJet(area_4vectors[parent2]; cluster_hist_index = i)
            elseif iszero(areas[parent2])
                area_4vectors[i] = PseudoJet(area_4vectors[parent1]; cluster_hist_index = i)
            else
                area_4vectors[i] = recombine(area_4vectors[parent1],
                                             area_4vectors[parent2];
                                             cluster_hist_index = i)
            end
        else
            error("Invalid cluster history at index $i: parent2 = $parent2") # COV_EXCL_LINE
        end
    end

    return JetAreaData(ghost_area,
                       n_ghosts,
                       n_hard_particles,
                       areas,
                       area_4vectors,
                       is_pure_ghost,
                       max_ghost_pt2,
                       has_dangerous_particles)
end

"""
    jet_reconstruct_area(particles; algorithm, R = 1.0,
                         area_spec = GhostedAreaSpec(), rng, strategy = RecoStrategy.Best,
                         p = nothing, recombine = addjets_escheme,
                         preprocess = preprocess_escheme)

Cluster particles with explicit ghosts and compute active jet areas.

Supports pp algorithms (`AntiKt`, `CA`, `Kt`, and `GenKt`, which requires `p`)
with all reconstruction strategies. The radius must be finite and positive.
For extreme `GenKt` powers, adjust `mean_ghost_pt` so ghost distance scales
neither overflow nor underflow; unrepresentable scales are rejected.
`rng` defaults to the task-local RNG, or a local seeded `MersenneTwister` when
`area_spec.seed` is set. Input particles are copied and assigned history indices.

Returns a [`ClusterSequenceArea`](@ref). Its `inclusive_jets` default to
`PseudoJet` and exclude pure ghost jets. Keep that type for area queries;
conversion to Lorentz vectors discards the history index.

This is a single explicit-ghost active-area measurement: `repeat` must be one,
and [`area_error`](@ref) is zero. Passive/Voronoi areas and pile-up subtraction
are not implemented. Check [`has_dangerous_particles`](@ref) if using unusually
soft physical particles or larger ghost momenta.
"""
function jet_reconstruct_area(particles::AbstractVector;
                              algorithm::JetAlgorithm.Algorithm,
                              R = 1.0,
                              p = nothing,
                              area_spec::GhostedAreaSpec = GhostedAreaSpec(),
                              rng::Random.AbstractRNG = _rng_from_area_spec(area_spec),
                              strategy::RecoStrategy.Strategy = RecoStrategy.Best,
                              recombine = addjets_escheme,
                              preprocess = preprocess_escheme,)
    if !is_pp(algorithm)
        throw(ArgumentError("jet_reconstruct_area currently supports pp algorithms only"))
    end
    radius = Float64(R)
    isfinite(radius) && radius > 0 && isfinite(radius * radius) && radius * radius > 0 ||
        throw(ArgumentError("R must have a finite, positive square"))
    power = get_algorithm_power(; algorithm, p)
    isfinite(power) || throw(ArgumentError("algorithm power must be finite"))

    hard_particles = Vector{PseudoJet}()
    sizehint!(hard_particles, length(particles))

    for (i, particle) in enumerate(particles)
        hard_particle = isnothing(preprocess) ?
                        PseudoJet(particle; cluster_hist_index = i) :
                        preprocess(particle, PseudoJet; cluster_hist_index = i)

        push!(hard_particles, hard_particle)
    end

    ghosts, actual_ghost_area = _make_area_ghosts(area_spec, rng)
    all(ghosts) do ghost
        distance_scale = pt2(ghost)^power * radius^2
        isfinite(distance_scale) && distance_scale > 0
    end ||
        throw(ArgumentError("ghost distance scales overflow or underflow; adjust mean_ghost_pt for this power and radius"))

    n_hard_particles = length(hard_particles)

    for i in eachindex(ghosts)
        ghosts[i] = PseudoJet(ghosts[i]; cluster_hist_index = n_hard_particles + i)
    end

    particles_with_ghosts = vcat(hard_particles, ghosts)

    cs = jet_reconstruct(particles_with_ghosts;
                         algorithm = algorithm,
                         R = radius,
                         p = p,
                         strategy = strategy,
                         recombine = recombine,
                         preprocess = nothing,)

    area_data = _postprocess_area(cs,
                                  n_hard_particles,
                                  actual_ghost_area;
                                  recombine = recombine,)

    return ClusterSequenceArea(cs, area_data)
end

"""Return the underlying `ClusterSequence`, including explicit ghosts."""
cluster_sequence(csa::ClusterSequenceArea) = csa.clusterseq

"""Return whether the area sequence contains explicit ghosts (always `true`)."""
has_explicit_ghosts(::ClusterSequenceArea) = true

function _area_history_index(jet::PseudoJet, csa::ClusterSequenceArea)
    index = cluster_hist_index(jet)
    cs = csa.clusterseq
    1 <= index <= length(cs.history) ||
        throw(ArgumentError("jet has no valid history index in this area sequence"))
    jet_index = cs.history[index].jetp_index
    1 <= jet_index <= length(cs.jets) && isequal(cs.jets[jet_index], jet) ||
        throw(ArgumentError("jet does not belong to this area sequence"))
    return index
end

"""Return the scalar active area of a `PseudoJet` belonging to `csa`."""
function area(jet::PseudoJet, csa::ClusterSequenceArea)
    return csa.area_data.areas[_area_history_index(jet, csa)]
end

"""Return the ghost area four-vector, propagated with the jet recombination scheme."""
function area_4vector(jet::PseudoJet, csa::ClusterSequenceArea)
    return csa.area_data.area_4vectors[_area_history_index(jet, csa)]
end

"""Return zero for a single explicit-ghost area measurement; no repeat error is estimated."""
function area_error(jet::PseudoJet, csa::ClusterSequenceArea)
    _area_history_index(jet, csa)
    return 0.0
end

"""Return whether a jet consists entirely of ghosts, independently of its momentum."""
function is_pure_ghost(jet::PseudoJet, csa::ClusterSequenceArea)
    return csa.area_data.is_pure_ghost[_area_history_index(jet, csa)]
end

"""Return the total scalar area of the generated ghost grid, including empty regions."""
function total_area(csa::ClusterSequenceArea)
    return csa.area_data.n_ghosts * csa.area_data.ghost_area
end

"""Return the largest squared transverse momentum of an initial ghost."""
function max_ghost_pt2(csa::ClusterSequenceArea)
    return csa.area_data.max_ghost_pt2
end

"""
    has_dangerous_particles(csa::ClusterSequenceArea)

Return whether any physical input has `eps(Float64)^2 * pt2 ≤ max_ghost_pt2`.
This is the FastJet softness diagnostic: ghosts may affect such particles'
clustering. Reduce ghost momenta or inspect the soft inputs when it is true.
"""
function has_dangerous_particles(csa::ClusterSequenceArea)
    return csa.area_data.has_dangerous_particles
end

"""Return the number of physical input particles, excluding added ghosts."""
function n_hard_particles(csa::ClusterSequenceArea)
    return csa.area_data.n_hard_particles
end

"""
    inclusive_jets(csa::ClusterSequenceArea, ::Type{T} = PseudoJet;
                   ptmin = 0.0, include_pure_ghosts = false)

Select inclusive jets, excluding pure ghost jets unless requested. The default
`PseudoJet` results retain the history needed by area queries. `LorentzVector`
and `LorentzVectorCyl` conversions are also supported but lose that history.
"""
function inclusive_jets(csa::ClusterSequenceArea,
                        ::Type{T} = PseudoJet;
                        ptmin = 0.0,
                        include_pure_ghosts::Bool = false,) where {T}
    T == PseudoJet || T <: LorentzVector || T <: LorentzVectorCyl ||
        throw(ArgumentError("unsupported return type $T for inclusive area jets"))
    jets = inclusive_jets(csa.clusterseq, PseudoJet; ptmin = ptmin)

    if !include_pure_ghosts
        filter!(jet -> !is_pure_ghost(jet, csa), jets)
    end

    T == PseudoJet && return jets

    if T <: LorentzVectorCyl
        return T[T(pt(jet), eta(jet), phi(jet), mass(jet)) for jet in jets]
    end
    return T[T(energy(jet), px(jet), py(jet), pz(jet)) for jet in jets]
end
