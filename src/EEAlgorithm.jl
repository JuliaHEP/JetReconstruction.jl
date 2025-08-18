# Use these constants whenever we need to set a large value for the distance or metric distance
const large_distance = 16.0 # = 4^2
const large_dij = 1.0e6

"""
    angular_distance(eereco, i, j) -> Float64

Calculate the angular distance between two jets `i` and `j` using the formula
``1 - cos(θ_{ij})``.

# Arguments
- `eereco`: The array of `EERecoJet` objects.
- `i`: The first jet.
- `j`: The second jet.

# Returns
- `Float64`: The angular distance between `i` and `j`, which is ``1 -
  cos\theta``.
"""
@inline function angular_distance(eereco, i, j)
    @inbounds @muladd 1.0 - eereco[i].nx * eereco[j].nx - eereco[i].ny * eereco[j].ny -
                      eereco[i].nz * eereco[j].nz
end

"""
    dij_dist(eereco, i, j, dij_factor, algorithm::JetAlgorithm.Algorithm, R = 4.0)

Calculate the dij distance between two e⁺e⁻ jets. This is the public entry point.
Internally, this forwards to a Val-based method for the given algorithm, which
allows the compiler to specialize away branches when `algorithm` is a constant.
"""
@inline function dij_dist(eereco, i, j, dij_factor, algorithm::JetAlgorithm.Algorithm,
                          R = 4.0)
    return dij_dist(eereco, i, j, dij_factor, Val(algorithm), R)
end

"""
    dij_dist(eereco, i, j, dij_factor, ::Val{JetAlgorithm.Durham}, R = 4.0)
    dij_dist(eereco, i, j, dij_factor, ::Val{JetAlgorithm.EEKt},   R = 4.0)

Durham/EEKt dij distance:
min(E_i^{2p}, E_j^{2p}) * dij_factor * (angular NN metric stored in nndist).
For EEKt, dij_factor encodes the R-dependent normalization.
"""
@inline function dij_dist(eereco, i, j, dij_factor, ::Val{JetAlgorithm.Durham}, R = 4.0)
    # Calculate the dij distance for jet i from jet j
    j == 0 && return large_dij
    @inbounds min(eereco[i].E2p, eereco[j].E2p) * dij_factor * eereco[i].nndist
end

@inline function dij_dist(eereco, i, j, dij_factor, ::Val{JetAlgorithm.EEKt}, R = 4.0)
    # Calculate the dij distance for jet i from jet j
    j == 0 && return large_dij
    @inbounds min(eereco[i].E2p, eereco[j].E2p) * dij_factor * eereco[i].nndist
end

# Fallback if a non-Algorithm token is passed
@inline function dij_dist(eereco, i, j, dij_factor, algorithm, R = 4.0)
    throw(ArgumentError("Algorithm $algorithm not supported for dij_dist"))
end

"""
    get_angular_nearest_neighbours!(eereco, ::Val{JetAlgorithm.Durham}, dij_factor, p, γ=1.0, R=4.0)

Compute initial nearest-neighbour distances for the Durham algorithm and fill
`eereco.nndist`, `eereco.nni` and `eereco.dijdist` accordingly. This variant
is specialised for the Durham metric and operates on the StructArray layout.
"""
@inline function get_angular_nearest_neighbours!(eereco, ::Val{JetAlgorithm.Durham},
                                                 dij_factor, p, γ = 1.0, R = 4.0)
    # Get the initial nearest neighbours for each jet
    N = length(eereco)
    # this_dist_vector = Vector{Float64}(undef, N)
    # Nearest neighbour geometric distance
    @inbounds for i in 1:N
        # TODO: Replace the 'j' loop with a vectorised operation over the appropriate array elements?
        # this_dist_vector .= 1.0 .- eereco.nx[i:N] .* eereco[i + 1:end].nx .-
        #     eereco[i].ny .* eereco[i + 1:end].ny .- eereco[i].nz .* eereco[i + 1:end].nz
        # The problem here will be avoiding allocations for the array outputs, which would easily
        # kill performance
        @inbounds for j in (i + 1):N
            this_nndist = angular_distance(eereco, i, j)

            # Using these ternary operators is faster than the if-else block
            better_nndist_i = this_nndist < eereco[i].nndist
            eereco.nndist[i] = better_nndist_i ? this_nndist : eereco.nndist[i]
            eereco.nni[i] = better_nndist_i ? j : eereco.nni[i]
            better_nndist_j = this_nndist < eereco[j].nndist
            eereco.nndist[j] = better_nndist_j ? this_nndist : eereco.nndist[j]
            eereco.nni[j] = better_nndist_j ? i : eereco.nni[j]
        end
    end
    # Nearest neighbour dij distance
    for i in 1:N
        eereco.dijdist[i] = dij_dist(eereco, i, eereco[i].nni, dij_factor,
                                     Val(JetAlgorithm.Durham), R)
    end
end

"""
    get_angular_nearest_neighbours!(eereco, ::Val{JetAlgorithm.EEKt}, dij_factor, p, γ=1.0, R=4.0)

Compute initial nearest-neighbour distances for the EEKt algorithm and fill
`eereco.nndist`, `eereco.nni` and `eereco.dijdist`. Performs beam-distance
checks appropriate for EEKt.
"""
@inline function get_angular_nearest_neighbours!(eereco, ::Val{JetAlgorithm.EEKt},
                                                 dij_factor, p, γ = 1.0, R = 4.0)
    # Get the initial nearest neighbours for each jet
    N = length(eereco)
    # this_dist_vector = Vector{Float64}(undef, N)
    # Nearest neighbour geometric distance
    @inbounds for i in 1:N
        # TODO: Replace the 'j' loop with a vectorised operation over the appropriate array elements?
        # this_dist_vector .= 1.0 .- eereco.nx[i:N] .* eereco[i + 1:end].nx .-
        #     eereco[i].ny .* eereco[i + 1:end].ny .- eereco[i].nz .* eereco[i + 1:end].nz
        # The problem here will be avoiding allocations for the array outputs, which would easily
        # kill performance
        @inbounds for j in (i + 1):N
            this_nndist = angular_distance(eereco, i, j)

            # Using these ternary operators is faster than the if-else block
            better_nndist_i = this_nndist < eereco[i].nndist
            eereco.nndist[i] = better_nndist_i ? this_nndist : eereco.nndist[i]
            eereco.nni[i] = better_nndist_i ? j : eereco.nni[i]
            better_nndist_j = this_nndist < eereco[j].nndist
            eereco.nndist[j] = better_nndist_j ? this_nndist : eereco.nndist[j]
            eereco.nni[j] = better_nndist_j ? i : eereco.nni[j]
        end
    end
    # Nearest neighbour dij distance
    for i in 1:N
        eereco.dijdist[i] = dij_dist(eereco, i, eereco[i].nni, dij_factor,
                                     Val(JetAlgorithm.EEKt), R)
    end
    # Beam-distance check for EEKt (this variant is EEKt-specific)
    @inbounds for i in 1:N
        beam_closer = eereco[i].E2p < eereco[i].dijdist
        eereco.dijdist[i] = beam_closer ? eereco[i].E2p : eereco.dijdist[i]
        eereco.nni[i] = beam_closer ? 0 : eereco.nni[i]
    end
end

"""
    get_angular_nearest_neighbours!(eereco, algorithm::JetAlgorithm.Algorithm, dij_factor, p, γ=1.0, R=4.0)

Forwarding wrapper that dispatches to the `Val{...}` specialised
`get_angular_nearest_neighbours!` implementation to avoid branches.
"""
@inline function get_angular_nearest_neighbours!(eereco,
                                                 algorithm::JetAlgorithm.Algorithm,
                                                 dij_factor, p, γ = 1.0, R = 4.0)
    return get_angular_nearest_neighbours!(eereco, Val(algorithm), dij_factor, p, γ, R)
end

"""
    update_nn_no_cross!(eereco, i, N, algorithm::JetAlgorithm.Algorithm, dij_factor, β=1.0, γ=1.0, R=4.0)

Update nearest-neighbour information for slot `i` (no-cross variant) by
forwarding to algorithm-specific `Val{...}` specialisations.
"""
@inline function update_nn_no_cross!(eereco, i, N, algorithm::JetAlgorithm.Algorithm,
                                     dij_factor, β = 1.0, γ = 1.0, R = 4.0)
    # Forward to Val-specialized implementations to avoid runtime branches
    return update_nn_no_cross!(eereco, i, N, Val(algorithm), dij_factor, β, γ, R)
end

# Val-specialized no-cross update
"""
    update_nn_no_cross!(eereco, i, N, ::Val{JetAlgorithm.Durham}, dij_factor, _β=1.0, _γ=1.0, R=4.0)

Durham-specialised update of the nearest neighbour for slot `i` (no-cross
variant). Updates `nndist`, `nni` and the inline `dijdist` for slot `i` using
direction-cosine arrays for performance.
"""
@inline function update_nn_no_cross!(eereco, i, N, ::Val{JetAlgorithm.Durham}, dij_factor,
                                     _β = 1.0, _γ = 1.0, R = 4.0)
    eereco.nndist[i] = large_distance
    eereco.nni[i] = i
    @inbounds for j in 1:N
        if j != i
            this_nndist = angular_distance(eereco, i, j)
            better_nndist_i = this_nndist < eereco[i].nndist
            eereco.nndist[i] = better_nndist_i ? this_nndist : eereco.nndist[i]
            eereco.nni[i] = better_nndist_i ? j : eereco.nni[i]
        end
    end
    eereco.dijdist[i] = dij_dist(eereco, i, eereco[i].nni, dij_factor,
                                 Val(JetAlgorithm.Durham), R)
end

"""
    update_nn_no_cross!(eereco, i, N, ::Val{JetAlgorithm.EEKt}, dij_factor, _β=1.0, _γ=1.0, R=4.0)

EEKt-specialised update of the nearest neighbour for slot `i` (no-cross
variant). Also applies the EEKt beam-distance check and updates `nni` when a
beam is closer than the dij distance.
"""
@inline function update_nn_no_cross!(eereco, i, N, ::Val{JetAlgorithm.EEKt}, dij_factor,
                                     _β = 1.0, _γ = 1.0, R = 4.0)
    # Update the nearest neighbour for jet i, w.r.t. all other active jets
    # also doing the cross check for the other jet
    eereco.nndist[i] = large_distance
    eereco.nni[i] = i
    @inbounds for j in 1:N
        if j != i
            this_nndist = angular_distance(eereco, i, j)
            better_nndist_i = this_nndist < eereco[i].nndist
            eereco.nndist[i] = better_nndist_i ? this_nndist : eereco.nndist[i]
            eereco.nni[i] = better_nndist_i ? j : eereco.nni[i]
            if this_nndist < eereco[j].nndist
                eereco.nndist[j] = this_nndist
                eereco.nni[j] = i
                # j will not be revisited, so update metric distance here
                eereco.dijdist[j] = dij_dist(eereco, j, i, dij_factor,
                                             Val(JetAlgorithm.EEKt), R)
                # EEKt-specific beam check is handled in the EEKt method only
            end
        end
    end
    eereco.dijdist[i] = dij_dist(eereco, i, eereco[i].nni, dij_factor,
                                 Val(JetAlgorithm.EEKt), R)

    # Need to check beam for EEKt
    beam_close = eereco[i].E2p < eereco[i].dijdist
    eereco.dijdist[i] = beam_close ? eereco[i].E2p : eereco.dijdist[i]
    eereco.nni[i] = beam_close ? 0 : eereco.nni[i]
end

"""
    update_nn_cross!(eereco, i, N, algorithm::JetAlgorithm.Algorithm, dij_factor, β=1.0, γ=1.0, R=4.0)

Forwarding wrapper for the cross-update variant. Dispatches to the
`Val{...}`-specialised `update_nn_cross!` to remove runtime branches.
"""
@inline function update_nn_cross!(eereco, i, N, algorithm::JetAlgorithm.Algorithm,
                                  dij_factor, β = 1.0, γ = 1.0, R = 4.0)
    # Forward to Val-specialized implementations to avoid runtime branches
    return update_nn_cross!(eereco, i, N, Val(algorithm), dij_factor, β, γ, R)
end

"""
    update_nn_cross!(eereco, i, N, ::Val{JetAlgorithm.Durham}, dij_factor, _β=1.0, _γ=1.0, R=4.0)

Durham-specialised cross-update: updates nearest-neighbour information for
slot `i` and propagates any changes to other slots when `i` becomes their new
nearest neighbour. Computes inline dij updates to avoid function-call overhead.
"""
@inline function update_nn_cross!(eereco, i, N, ::Val{JetAlgorithm.Durham}, dij_factor,
                                  _β = 1.0, _γ = 1.0, R = 4.0)
    # Update the nearest neighbour for jet i, w.r.t. all other active jets
    # also doing the cross check for the other jet
    eereco.nndist[i] = large_distance
    eereco.nni[i] = i
    @inbounds for j in 1:N
        if j != i
            this_nndist = angular_distance(eereco, i, j)
            better_nndist_i = this_nndist < eereco[i].nndist
            eereco.nndist[i] = better_nndist_i ? this_nndist : eereco.nndist[i]
            eereco.nni[i] = better_nndist_i ? j : eereco.nni[i]
            if this_nndist < eereco[j].nndist
                eereco.nndist[j] = this_nndist
                eereco.nni[j] = i
                # j will not be revisited, so update metric distance here
                eereco.dijdist[j] = dij_dist(eereco, j, i, dij_factor,
                                             Val(JetAlgorithm.Durham), R)
                # EEKt-specific beam check is handled in the EEKt method only
            end
        end
    end
    eereco.dijdist[i] = dij_dist(eereco, i, eereco[i].nni, dij_factor,
                                 Val(JetAlgorithm.Durham), R)
end

"""
    update_nn_cross!(eereco, i, N, ::Val{JetAlgorithm.EEKt}, dij_factor, _β=1.0, _γ=1.0, R=4.0)

EEKt-specialised cross-update. Updates `nndist`, `nni` and `dijdist` for
both the updated slot and any affected neighbours, performing beam checks
where required.
"""
@inline function update_nn_cross!(eereco, i, N, ::Val{JetAlgorithm.EEKt}, dij_factor,
                                  _β = 1.0, _γ = 1.0, R = 4.0)
    # Update the nearest neighbour for jet i, w.r.t. all other active jets
    # also doing the cross check for the other jet
    eereco.nndist[i] = large_distance
    eereco.nni[i] = i
    @inbounds for j in 1:N
        if j != i
            this_nndist = angular_distance(eereco, i, j)
            better_nndist_i = this_nndist < eereco[i].nndist
            eereco.nndist[i] = better_nndist_i ? this_nndist : eereco.nndist[i]
            eereco.nni[i] = better_nndist_i ? j : eereco.nni[i]
            if this_nndist < eereco[j].nndist
                eereco.nndist[j] = this_nndist
                eereco.nni[j] = i
                # j will not be revisited, so update metric distance here
                eereco.dijdist[j] = dij_dist(eereco, j, i, dij_factor,
                                             Val(JetAlgorithm.EEKt), R)
                if eereco[j].E2p < eereco[j].dijdist
                    eereco.dijdist[j] = eereco[j].E2p
                    eereco.nni[j] = 0
                end
            end
        end
    end
    eereco.dijdist[i] = dij_dist(eereco, i, eereco[i].nni, dij_factor,
                                 Val(JetAlgorithm.EEKt), R)

    # Check beam for EEKt
    beam_close = eereco[i].E2p < eereco[i].dijdist
    eereco.dijdist[i] = beam_close ? eereco[i].E2p : eereco.dijdist[i]
    eereco.nni[i] = beam_close ? 0 : eereco.nni[i]
end

"""
    ee_check_consistency(clusterseq, eereco, N)

Run an internal consistency check over the reconstruction state. Emits
`@error` logs for invalid nearest-neighbour indices or malformed history
entries. Intended for debugging; inexpensive checks only.
"""
function ee_check_consistency(clusterseq, eereco, N)
    # Check the consistency of the reconstruction state
    for i in 1:N
        if eereco[i].nni > N
            @error "Jet $i has invalid nearest neighbour $(eereco[i].nni)"
        end
        for i in 1:N
            jet = clusterseq.jets[eereco[i].index]
            jet_hist = clusterseq.history[jet._cluster_hist_index]
            if jet_hist.child != Invalid
                @error "Jet $i has invalid child $(jet_hist.child)"
            end
        end
    end
    @debug "Consistency check passed"
end

Base.@propagate_inbounds @inline function fill_reco_array!(eereco, particles, R2, p)
    @inbounds for i in eachindex(particles)
        eereco.index[i] = i
        eereco.nni[i] = 0
        eereco.nndist[i] = R2
        # eereco.dijdist[i] = UNDEF # Does not need to be initialised
        eereco.nx[i] = nx(particles[i])
        eereco.ny[i] = ny(particles[i])
        eereco.nz[i] = nz(particles[i])
        eereco.E2p[i] = energy(particles[i])^(2p)
        # No precomputed beam factor; compute on demand to preserve previous behaviour
    end
end

Base.@propagate_inbounds @inline function insert_new_jet!(eereco, i, newjet_k, R2,
                                                          merged_jet, p)
    @inbounds begin
        eereco.index[i] = newjet_k
        eereco.nni[i] = 0
        eereco.nndist[i] = R2
        eereco.nx[i] = nx(merged_jet)
        eereco.ny[i] = ny(merged_jet)
        eereco.nz[i] = nz(merged_jet)
        E = energy(merged_jet)
        if p isa Int
            if p == 1
                eereco.E2p[i] = E * E
            else
                E2 = E * E
                eereco.E2p[i] = E2^p
            end
        else
            eereco.E2p[i] = E^(2p)
        end
    end
end

"""
    copy_to_slot!(eereco, i, j)

Copy the contents of slot `i` in the `eereco` array to slot `j`.
"""
Base.@propagate_inbounds @inline function copy_to_slot!(eereco, i, j)
    @inbounds begin
        eereco.index[j] = eereco.index[i]
        eereco.nni[j] = eereco.nni[i]
        eereco.nndist[j] = eereco.nndist[i]
        eereco.dijdist[j] = eereco.dijdist[i]
        eereco.nx[j] = eereco.nx[i]
        eereco.ny[j] = eereco.ny[i]
        eereco.nz[j] = eereco.nz[i]
        eereco.E2p[j] = eereco.E2p[i]
    end
end

"""
    ee_genkt_algorithm(particles::AbstractVector{T}; algorithm::JetAlgorithm.Algorithm,
                       p::Union{Real, Nothing} = nothing, R = 4.0, recombine = addjets,
                       preprocess = nothing, γ::Real = 1.0) where {T}

Run an e+e- reconstruction algorithm on a set of initial particles.

# Arguments
- `particles::AbstractVector{T}`: A vector of particles to be clustered.
- `algorithm::JetAlgorithm.Algorithm`: The jet algorithm to use.
- `p::Union{Real, Nothing} = nothing`: The power parameter for the algorithm.
  Must be specified for EEKt algorithm. 
  For Valencia algorithm, this corresponds to the β parameter.
  Other algorithms will ignore this value.
- `R = 4.0`: The jet radius parameter. Not required / ignored for the Durham
  algorithm.
- `recombine`: The recombination scheme to use.
- `preprocess`: Preprocessing function for input particles.
- `γ::Real = 1.0`: The angular exponent parameter for Valencia algorithm. Ignored for other algorithms.

# Returns
- The result of the jet clustering as a `ClusterSequence` object.

# Notes
This is the public interface to the e+e- jet clustering algorithm. The function
will check for consistency between the algorithm and the power parameter as
needed. It will then prepare the internal EDM particles for the clustering
itself, and call the actual reconstruction method `_ee_genkt_algorithm!`.

If the algorithm is Durham, `R` is nominally set to 4.
If the algorithm is EEkt, power `p` must be specified.
If the algorithm is Valencia, both `p` (β) and `γ` should be specified.
"""
function ee_genkt_algorithm(particles::AbstractVector{T}; algorithm::JetAlgorithm.Algorithm,
                            p::Union{Real, Nothing} = nothing, R = 4.0, recombine = addjets,
                            preprocess = nothing, γ::Real = 1.0,
                            β::Union{Real, Nothing} = nothing) where {T}
    # (β override for Valencia handled inside _ee_genkt_algorithm_valencia)
    if algorithm == JetAlgorithm.Valencia
        # Apply β override if provided, then obtain validated power
        local_p = β === nothing ? p : β
        local_p = get_algorithm_power(p = local_p, algorithm = algorithm)
        return _ee_genkt_algorithm_valencia(particles = begin
                                                if isnothing(preprocess)
                                                    if T == EEJet
                                                        recombination_particles = copy(particles)
                                                        sizehint!(recombination_particles,
                                                                  length(particles) * 2)
                                                        recombination_particles
                                                    else
                                                        recombination_particles = EEJet[]
                                                        sizehint!(recombination_particles,
                                                                  length(particles) * 2)
                                                        for (i, particle) in enumerate(particles)
                                                            push!(recombination_particles,
                                                                  EEJet(particle;
                                                                        cluster_hist_index = i))
                                                        end
                                                        recombination_particles
                                                    end
                                                else
                                                    recombination_particles = EEJet[]
                                                    sizehint!(recombination_particles,
                                                              length(particles) * 2)
                                                    for (i, particle) in enumerate(particles)
                                                        push!(recombination_particles,
                                                              preprocess(particle, EEJet;
                                                                         cluster_hist_index = i))
                                                    end
                                                    recombination_particles
                                                end
                                            end,
                                            algorithm = algorithm, p = local_p, R = R,
                                            recombine = recombine, γ = γ, beta = β)
    end

    # Check for consistency algorithm power (non-Valencia path)
    p = get_algorithm_power(p = p, algorithm = algorithm)

    # Cast p to Int where possible for algorithms that benefit (Durham/EEKt)
    if algorithm == JetAlgorithm.Durham || algorithm == JetAlgorithm.EEKt
        p = (round(p) == p) ? Int(p) : p
    end

    # For the Durham algorithm, p=1 and R is not used, but nominally set to 4
    if algorithm == JetAlgorithm.Durham
        R = 4.0
    end

    if isnothing(preprocess)
        if T == EEJet
            # If we don't have a preprocessor, we just need to copy to our own
            # EEJet objects
            recombination_particles = copy(particles)
            sizehint!(recombination_particles, length(particles) * 2)
        else
            # We assume a constructor for EEJet that can ingest the appropriate
            # type of particle
            recombination_particles = EEJet[]
            sizehint!(recombination_particles, length(particles) * 2)
            for (i, particle) in enumerate(particles)
                push!(recombination_particles, EEJet(particle; cluster_hist_index = i))
            end
        end
    else
        # We have a preprocessor function that we need to call to modify the
        # input particles
        recombination_particles = EEJet[]
        sizehint!(recombination_particles, length(particles) * 2)
        for (i, particle) in enumerate(particles)
            push!(recombination_particles,
                  preprocess(particle, EEJet; cluster_hist_index = i))
        end
    end

    # Now call the actual reconstruction method.
    if algorithm == JetAlgorithm.Valencia
        return _ee_genkt_algorithm_valencia(particles = recombination_particles, p = p,
                                            R = R,
                                            algorithm = algorithm,
                                            recombine = recombine, γ = γ)
    elseif algorithm == JetAlgorithm.Durham
        return _ee_genkt_algorithm_durham(particles = recombination_particles, p = p, R = R,
                                          algorithm = algorithm,
                                          recombine = recombine, γ = γ)
    else
        return _ee_genkt_algorithm(particles = recombination_particles, p = p, R = R,
                                   algorithm = algorithm,
                                   recombine = recombine, γ = γ)
    end
end

################################################################################
# Durham-specialised implementation
################################################################################
function _ee_genkt_algorithm_durham(; particles::AbstractVector{EEJet},
                                    algorithm::JetAlgorithm.Algorithm, p::Real, R = 4.0,
                                    recombine = addjets, γ::Real = 1.0,
                                    beta::Union{Real, Nothing} = nothing)
    # Bounds
    N::Int = length(particles)

    R2 = R^2

    # Durham dij factor
    dij_factor = 2.0

    # Prepare SoA
    eereco = StructArray{EERecoJet}(undef, N)
    fill_reco_array!(eereco, particles, R2, p)

    # Setup history
    history, Qtot = initial_history(particles)
    clusterseq = ClusterSequence(algorithm, p, R, RecoStrategy.N2Plain, particles, history,
                                 Qtot)

    # Initial nearest neighbours (Durham-specialised)
    get_angular_nearest_neighbours!(eereco, Val(JetAlgorithm.Durham), dij_factor, p, γ, R)

    # Alias StructArray fields to local vectors to avoid repeated getindex
    index = eereco.index
    nni = eereco.nni
    nndist = eereco.nndist
    dijdist = eereco.dijdist
    nx = eereco.nx
    ny = eereco.ny
    nz = eereco.nz
    E2p = eereco.E2p

    # Main loop
    iter = 0
    @inbounds while N != 0
        iter += 1

        dij_min, ijetA = fast_findmin(dijdist, N)
        ijetB = nni[ijetA]

        if ijetB == 0
            ijetB = ijetA
            add_step_to_history!(clusterseq,
                                 clusterseq.jets[index[ijetA]]._cluster_hist_index,
                                 BeamJet, Invalid, dij_min)
        elseif N == 1
            ijetB = ijetA
            add_step_to_history!(clusterseq,
                                 clusterseq.jets[eereco[ijetA].index]._cluster_hist_index,
                                 BeamJet, Invalid, dij_min)
        else
            if ijetB < ijetA
                ijetA, ijetB = ijetB, ijetA
            end
            jetA = clusterseq.jets[index[ijetA]]
            jetB = clusterseq.jets[index[ijetB]]
            merged_jet = recombine(jetA, jetB;
                                   cluster_hist_index = length(clusterseq.history) + 1)
            push!(clusterseq.jets, merged_jet)
            newjet_k = length(clusterseq.jets)
            add_step_to_history!(clusterseq,
                                 minmax(cluster_hist_index(jetA),
                                        cluster_hist_index(jetB))...,
                                 newjet_k, dij_min)
            insert_new_jet!(eereco, ijetA, newjet_k, R2, merged_jet, p)
        end

        if ijetB != N
            copy_to_slot!(eereco, N, ijetB)
        end

        N -= 1

        # Update nearest neighbours step using array-based helpers specialised for Durham
        for i in 1:N
            if (ijetB != N + 1) && (nni[i] == N + 1)
                nni[i] = ijetB
            else
                if (nni[i] == ijetA) || (nni[i] == ijetB) || (nni[i] > N)
                    # use SoA-based updater
                    update_nn_no_cross!(eereco, i, N, Val(JetAlgorithm.Durham), dij_factor,
                                        p, γ, R)
                end
            end
        end

        if ijetA != ijetB
            update_nn_cross!(eereco, ijetA, N, Val(JetAlgorithm.Durham), dij_factor, p, γ,
                             R)
        end
    end

    clusterseq
end

################################################################################
# EEKt-specialised implementation
################################################################################

"""
    _ee_genkt_algorithm!(particles::AbstractVector{EEJet};
                        algorithm::JetAlgorithm.Algorithm, p::Real, R = 4.0,
                        recombine = addjets, γ::Real = 1.0)

This function is the internal implementation of the e+e- jet clustering
algorithm. It takes a vector of `EEJet` `particles` representing the input
particles and reconstructs jets based on the specified parameters.

Users of the package should use the `ee_genkt_algorithm` function as their
entry point to this jet reconstruction.

# Arguments
- `particles::AbstractVector{EEJet}`: A vector of `EEJet` particles used
  as input for jet reconstruction. This vector must supply the correct
  `cluster_hist_index` values and will be *mutated* as part of the returned
  `ClusterSequence`.
- `algorithm::JetAlgorithm.Algorithm`: The jet reconstruction algorithm to use.
- `p::Real`: The power to which the transverse momentum (`pt`) of each particle
  is raised.
- `R = 4.0`: The jet radius parameter.
- `recombine = addjets`: The recombination function used to merge two jets.

# Returns
- `clusterseq`: The resulting `ClusterSequence` object representing the
  reconstructed jets.
"""
function _ee_genkt_algorithm(; particles::AbstractVector{EEJet},
                             algorithm::JetAlgorithm.Algorithm, p::Real, R = 4.0,
                             recombine = addjets, γ::Real = 1.0,
                             beta::Union{Real, Nothing} = nothing)
    # Bounds
    N::Int = length(particles)

    # Forward Valencia directly to specialised implementation (tests call this internal API)
    if algorithm == JetAlgorithm.Valencia
        return _ee_genkt_algorithm_valencia(particles = particles, algorithm = algorithm,
                                            p = (beta === nothing ? p : beta), R = R,
                                            recombine = recombine, γ = γ, beta = beta)
    end

    R2 = R^2

    # Constant factor for the dij metric and the beam distance function
    if algorithm == JetAlgorithm.Durham
        dij_factor = 2.0
    elseif algorithm == JetAlgorithm.EEKt
        if R < π
            dij_factor = 1 / (1 - cos(R))
        else
            dij_factor = 1 / (3 + cos(R))
        end
    else
        throw(ArgumentError("Algorithm $algorithm not supported for e+e-"))
    end

    # For optimised reconstruction generate an SoA containing the necessary
    # jet information and populate it accordingly
    # We need N slots for this array
    eereco = StructArray{EERecoJet}(undef, N)

    fill_reco_array!(eereco, particles, R2, p)

    # Setup the initial history and get the total energy
    history, Qtot = initial_history(particles)

    clusterseq = ClusterSequence(algorithm, p, R, RecoStrategy.N2Plain, particles, history,
                                 Qtot)

    # Run over initial pairs of jets to find nearest neighbours
    get_angular_nearest_neighbours!(eereco, algorithm, dij_factor, p, γ, R)

    # Only for debugging purposes...
    # ee_check_consistency(clusterseq, clusterseq_index, N, nndist, nndij, nni, "Start")

    # Now we can start the main loop
    iter = 0
    while N != 0
        iter += 1

        dij_min, ijetA = fast_findmin(eereco.dijdist, N)
        ijetB = eereco[ijetA].nni

        # Now we check if there is a "beam" merge possibility
        if ijetB == 0
            # We have an EEKt beam merge
            ijetB = ijetA
            add_step_to_history!(clusterseq,
                                 clusterseq.jets[eereco[ijetA].index]._cluster_hist_index,
                                 BeamJet, Invalid, dij_min)
        elseif N == 1
            # We have a single jet left
            ijetB = ijetA
            add_step_to_history!(clusterseq,
                                 clusterseq.jets[eereco[ijetA].index]._cluster_hist_index,
                                 BeamJet, Invalid, dij_min)
        else
            # Jet-jet merge
            if ijetB < ijetA
                ijetA, ijetB = ijetB, ijetA
            end

            # Resolve the jet indexes to access the actual jets
            jetA = clusterseq.jets[eereco[ijetA].index]
            jetB = clusterseq.jets[eereco[ijetB].index]

            # Recombine jetA and jetB into the next jet
            merged_jet = recombine(jetA, jetB;
                                   cluster_hist_index = length(clusterseq.history) + 1)

            # Now add the jet to the sequence, and update the history
            push!(clusterseq.jets, merged_jet)
            newjet_k = length(clusterseq.jets)
            add_step_to_history!(clusterseq,
                                 minmax(cluster_hist_index(jetA),
                                        cluster_hist_index(jetB))...,
                                 newjet_k, dij_min)

            # Update the compact arrays, reusing the JetA slot
            insert_new_jet!(eereco, ijetA, newjet_k, R2, merged_jet, p)
        end

        # Squash step - copy the final jet's compact data into the jetB slot
        # unless we are at the end of the array, in which case do nothing
        if ijetB != N
            copy_to_slot!(eereco, N, ijetB)
        end

        # Now number of active jets is decreased by one
        N -= 1

        # Update nearest neighbours step
        for i in 1:N
            # First, if any jet's NN was the old index N, it's now ijetB
            if (ijetB != N + 1) && (eereco[i].nni == N + 1)
                eereco.nni[i] = ijetB
            else
                # Otherwise, if the jet had ijetA or ijetB as their NN, we need to update them
                # plus "belt and braces" check for an invalid NN (>N)
                if (eereco[i].nni == ijetA) || (eereco[i].nni == ijetB) ||
                   (eereco[i].nni > N)
                    update_nn_no_cross!(eereco, i, N, algorithm, dij_factor, p, γ, R)
                end
            end
        end

        # Finally, we need to update the nearest neighbours for the new jet, checking both ways
        # (But only if there was a new jet!)
        if ijetA != ijetB
            update_nn_cross!(eereco, ijetA, N, algorithm, dij_factor, p, γ, R)
        end

        # Only for debugging purposes...
        # ee_check_consistency(clusterseq, eereco, N)
    end

    # Return the final cluster sequence structure
    clusterseq
end
