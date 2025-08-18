################################################################################
# Valencia-specialised helpers and implementation
################################################################################

"""
    valencia_distance_inv(eereco, i, j, invR2) -> Float64

Calculate the Valencia dij metric (scaled by `invR2`) between slots `i` and
`j` in a StructArray `eereco`. Uses E2p and direction cosines from `eereco`.
"""
Base.@propagate_inbounds @inline function valencia_distance_inv(eereco, i, j, invR2)
    # Assume SoA layout
    nx = eereco.nx
    ny = eereco.ny
    nz = eereco.nz
    E2p = eereco.E2p
    angular_dist = angular_distance_arrays(nx, ny, nz, i, j)
    # Valencia dij : min(E_i^{2β}, E_j^{2β}) * 2 * (1 - cos θ) * invR2
    min(E2p[i], E2p[j]) * 2 * angular_dist * invR2
end

"""
    dij_dist(eereco, i, j, dij_factor, ::Val{JetAlgorithm.Valencia}, R)

Valencia-specialised `dij_dist` which computes the Valencia dij metric for
slots `i` and `j`. Returns a large sentinel distance for beam index `j==0`.
"""
@inline function dij_dist(eereco, i, j, dij_factor, ::Val{JetAlgorithm.Valencia}, R)
    j == 0 && return large_dij
    @inbounds valencia_distance(eereco, i, j, R)
end

"""
    valencia_distance(eereco, i, j, R) -> Float64

Compute the Valencia dij metric between `i` and `j` using explicit `R` and
delegating to `valencia_distance_inv` with the appropriate scaling.
"""
Base.@propagate_inbounds @inline function valencia_distance(eereco, i, j, R)
    return valencia_distance_inv(eereco, i, j, inv(R * R))
end

# Array-based helpers
"""
    valencia_distance_inv_arrays(E2p, nx, ny, nz, i, j, invR2) -> Float64

Array-based variant of `valencia_distance_inv` that works directly on raw
vectors (useful for the precomputed helpers/fast paths).
"""
Base.@propagate_inbounds @inline function valencia_distance_inv_arrays(E2p, nx, ny, nz, i,
                                                                       j, invR2)
    angular_dist = angular_distance_arrays(nx, ny, nz, i, j)
    min(E2p[i], E2p[j]) * 2 * angular_dist * invR2
end

# Scaled variant
"""
    valencia_distance_inv_scaled_arrays(E2p_scaled, nx, ny, nz, i, j) -> Float64

Compute the Valencia distance using a pre-scaled E2p vector (`E2p_scaled`),
avoiding repeated multiplication by the R-dependent factor. Intended for
performance-sensitive precomputed loops.
"""
Base.@propagate_inbounds @inline function valencia_distance_inv_scaled_arrays(E2p_scaled,
                                                                              nx, ny, nz, i,
                                                                              j)
    angular_dist = angular_distance_arrays(nx, ny, nz, i, j)
    min(E2p_scaled[i], E2p_scaled[j]) * angular_dist
end

"""
    valencia_beam_distance(eereco, i, γ, β) -> Float64

Compute the Valencia beam-distance term for slot `i` given angular exponent
`γ` and (unused) `β`. Uses direction cosine `nz` and E2p value.
"""
Base.@propagate_inbounds @inline function valencia_beam_distance(eereco, i, γ, β)
    nzv = eereco.nz
    E2pv = eereco.E2p
    nz_i = nzv[i]
    sin2 = 1 - nz_i * nz_i
    E2p = E2pv[i]
    if γ == 1.0
        return E2p * sin2
    elseif γ == 2.0
        return E2p * (sin2 * sin2)
    else
        return E2p * sin2^γ
    end
end

"""
    valencia_beam_distance_arrays(E2p, nz, i, γ, β) -> Float64

Array-based variant of `valencia_beam_distance` that operates on raw `E2p`
and `nz` vectors.
"""
Base.@propagate_inbounds @inline function valencia_beam_distance_arrays(E2p, nz, i, γ, β)
    nz_i = nz[i]
    sin2 = 1 - nz_i * nz_i
    E2p_i = E2p[i]
    if γ == 1.0
        return E2p_i * sin2
    elseif γ == 2.0
        return E2p_i * (sin2 * sin2)
    else
        return E2p_i * sin2^γ
    end
end

"""
    get_angular_nearest_neighbours_valencia_precomputed!(eereco, E2p_scaled, beam_term, p, γ=1.0, R=4.0)

Initialize nearest-neighbour arrays for the Valencia algorithm using
precomputed scaled energy vector `E2p_scaled` and `beam_term`. This avoids
repeated per-pair scaling inside tight loops.
"""
Base.@propagate_inbounds @inline function get_angular_nearest_neighbours_valencia_precomputed!(eereco,
                                                                                               E2p_scaled::AbstractVector,
                                                                                               beam_term::AbstractVector,
                                                                                               p,
                                                                                               γ = 1.0,
                                                                                               R = 4.0)
    N = length(eereco)
    nx = eereco.nx
    ny = eereco.ny
    nz = eereco.nz
    nndist = eereco.nndist
    nni = eereco.nni
    @inbounds for i in 1:N
        nndist[i] = Inf
        nni[i] = i
    end
    @inbounds for i in 1:N
        local_nndist_i = nndist[i]
        local_nni_i = nni[i]
        @inbounds for j in (i + 1):N
            this_metric = valencia_distance_inv_scaled_arrays(E2p_scaled, nx, ny, nz, i, j)
            if this_metric < local_nndist_i
                local_nndist_i = this_metric
                local_nni_i = j
            end
            if this_metric < nndist[j]
                nndist[j] = this_metric
                nni[j] = i
            end
        end
        nndist[i] = local_nndist_i
        nni[i] = local_nni_i
    end
    @inbounds for i in 1:N
        eereco.dijdist[i] = valencia_distance_inv_scaled_arrays(E2p_scaled, nx, ny, nz, i,
                                                                nni[i])
    end
    @inbounds for i in 1:N
        beam_closer = beam_term[i] < eereco.dijdist[i]
        eereco.dijdist[i] = beam_closer ? beam_term[i] : eereco.dijdist[i]
        eereco.nni[i] = beam_closer ? 0 : eereco.nni[i]
    end
end

"""
    update_nn_no_cross_arrays_precomputed!(nndist, nni, nx, ny, nz, E2p_scaled, beam_term, dijdist, i, N, dij_factor, β=1.0, γ=1.0, R=4.0)

Precomputed Valencia no-cross nearest-neighbour update. Uses `E2p_scaled`
and `beam_term` arrays to compute distances without per-pair scaling.
"""
@inline function update_nn_no_cross_arrays_precomputed!(nndist::AbstractVector,
                                                        nni::AbstractVector,
                                                        nx::AbstractVector,
                                                        ny::AbstractVector,
                                                        nz::AbstractVector,
                                                        E2p_scaled::AbstractVector,
                                                        beam_term::AbstractVector,
                                                        dijdist::AbstractVector,
                                                        i::Integer, N::Integer,
                                                        dij_factor, β = 1.0, γ = 1.0,
                                                        R = 4.0)
    nndist_i = Inf
    nni_i = i
    @inbounds for j in 1:N
        if j != i
            this_metric = valencia_distance_inv_scaled_arrays(E2p_scaled, nx, ny, nz, i, j)
            if this_metric < nndist_i
                nndist_i = this_metric
                nni_i = j
            end
        end
    end
    dijdist_i = valencia_distance_inv_scaled_arrays(E2p_scaled, nx, ny, nz, i, nni_i)
    if beam_term[i] < dijdist_i
        dijdist_i = beam_term[i]
        nni_i = 0
    end
    nndist[i] = nndist_i
    nni[i] = nni_i
    dijdist[i] = dijdist_i
end

"""
    update_nn_cross_arrays_precomputed!(nndist, nni, nx, ny, nz, E2p_scaled, beam_term, dijdist, i, N, dij_factor, β=1.0, γ=1.0, R=4.0)

Precomputed Valencia cross-update variant: updates neighbour data for slot
`i` and any affected neighbours using `E2p_scaled` and `beam_term`.
"""
@inline function update_nn_cross_arrays_precomputed!(nndist::AbstractVector,
                                                     nni::AbstractVector,
                                                     nx::AbstractVector, ny::AbstractVector,
                                                     nz::AbstractVector,
                                                     E2p_scaled::AbstractVector,
                                                     beam_term::AbstractVector,
                                                     dijdist::AbstractVector,
                                                     i::Integer, N::Integer,
                                                     dij_factor, β = 1.0, γ = 1.0, R = 4.0)
    nndist_i = Inf
    nni_i = i
    @inbounds for j in 1:(i - 1)
        this_metric = valencia_distance_inv_scaled_arrays(E2p_scaled, nx, ny, nz, i, j)
        if this_metric < nndist_i
            nndist_i = this_metric
            nni_i = j
        end
        if this_metric < nndist[j]
            nndist[j] = this_metric
            nni[j] = i
            dijdist[j] = this_metric
            if beam_term[j] < dijdist[j]
                dijdist[j] = beam_term[j]
                nni[j] = 0
            end
        end
    end
    @inbounds for j in (i + 1):N
        this_metric = valencia_distance_inv_scaled_arrays(E2p_scaled, nx, ny, nz, i, j)
        if this_metric < nndist_i
            nndist_i = this_metric
            nni_i = j
        end
        if this_metric < nndist[j]
            nndist[j] = this_metric
            nni[j] = i
            dijdist[j] = this_metric
            if beam_term[j] < dijdist[j]
                dijdist[j] = beam_term[j]
                nni[j] = 0
            end
        end
    end
    # Finalize slot i
    dijdist_i = valencia_distance_inv_scaled_arrays(E2p_scaled, nx, ny, nz, i, nni_i)
    if beam_term[i] < dijdist_i
        dijdist_i = beam_term[i]
        nni_i = 0
    end
    nndist[i] = nndist_i
    nni[i] = nni_i
    dijdist[i] = dijdist_i
end

"""
    _ee_genkt_algorithm_valencia(; particles::AbstractVector{EEJet}, algorithm::JetAlgorithm.Algorithm, p::Real, R=4.0, recombine=addjets, γ::Real=1.0, beta::Union{Real, Nothing}=nothing)

Valencia-specialised implementation of the e+e- gen-kT clustering algorithm.
This implementation precomputes scaled energy and beam-term arrays to speed up
nearest-neighbour computations for the Valencia metric.
"""
function _ee_genkt_algorithm_valencia(; particles::AbstractVector{EEJet},
                                      algorithm::JetAlgorithm.Algorithm, p::Real, R = 4.0,
                                      recombine = addjets, γ::Real = 1.0,
                                      beta::Union{Real, Nothing} = nothing)
    N::Int = length(particles)
    R2 = R^2
    if algorithm == JetAlgorithm.Valencia && beta !== nothing
        p = beta
    end
    dij_factor = 1.0
    eereco = StructArray{EERecoJet}(undef, N)
    fill_reco_array!(eereco, particles, R2, p)
    history, Qtot = initial_history(particles)
    clusterseq = ClusterSequence(algorithm, p, R, RecoStrategy.N2Plain, particles, history,
                                 Qtot)
    indexv = eereco.index
    nni_v = eereco.nni
    nndist_v = eereco.nndist
    dijdist_v = eereco.dijdist
    nxv = eereco.nx
    nyv = eereco.ny
    nzv = eereco.nz
    E2pv = eereco.E2p
    invR2 = inv(R * R)
    factor = 2 * invR2
    E2p_scaled = similar(E2pv)
    beam_term = similar(E2pv)
    @inbounds for k in 1:N
        E2p_scaled[k] = E2pv[k] * factor
        nz_k = nzv[k]
        sin2 = 1.0 - nz_k * nz_k
        if γ == 1.0
            beam_term[k] = E2pv[k] * sin2
        elseif γ == 2.0
            beam_term[k] = E2pv[k] * (sin2 * sin2)
        else
            beam_term[k] = E2pv[k] * sin2^γ
        end
    end
    get_angular_nearest_neighbours_valencia_precomputed!(eereco, E2p_scaled, beam_term, p,
                                                         γ, R)
    iter = 0
    while N != 0
        iter += 1
        dij_min, ijetA = fast_findmin(dijdist_v, N)
        ijetB = nni_v[ijetA]
        if ijetB == 0
            ijetB = ijetA
            add_step_to_history!(clusterseq,
                                 clusterseq.jets[indexv[ijetA]]._cluster_hist_index,
                                 BeamJet, Invalid, dij_min)
        elseif N == 1
            ijetB = ijetA
            add_step_to_history!(clusterseq,
                                 clusterseq.jets[indexv[ijetA]]._cluster_hist_index,
                                 BeamJet, Invalid, dij_min)
        else
            if ijetB < ijetA
                ijetA, ijetB = ijetB, ijetA
            end
            jetA = clusterseq.jets[indexv[ijetA]]
            jetB = clusterseq.jets[indexv[ijetB]]
            merged_jet = recombine(jetA, jetB;
                                   cluster_hist_index = length(clusterseq.history) + 1)
            push!(clusterseq.jets, merged_jet)
            newjet_k = length(clusterseq.jets)
            add_step_to_history!(clusterseq,
                                 minmax(cluster_hist_index(jetA),
                                        cluster_hist_index(jetB))...,
                                 newjet_k, dij_min)
            indexv[ijetA] = newjet_k
            nni_v[ijetA] = 0
            nndist_v[ijetA] = R2
            nxv[ijetA] = nx(merged_jet)
            nyv[ijetA] = ny(merged_jet)
            nzv[ijetA] = nz(merged_jet)
            E = energy(merged_jet)
            if p isa Int
                if p == 1
                    E2pv[ijetA] = E * E
                else
                    E2 = E * E
                    E2pv[ijetA] = E2^p
                end
            else
                E2pv[ijetA] = E^(2p)
            end
            E2p_scaled[ijetA] = E2pv[ijetA] * factor
            nz_k = nzv[ijetA]
            sin2_k = 1.0 - nz_k * nz_k
            if γ == 1.0
                beam_term[ijetA] = E2pv[ijetA] * sin2_k
            elseif γ == 2.0
                beam_term[ijetA] = E2pv[ijetA] * (sin2_k * sin2_k)
            else
                beam_term[ijetA] = E2pv[ijetA] * sin2_k^γ
            end
        end
        if ijetB != N
            indexv[ijetB] = indexv[N]
            nni_v[ijetB] = nni_v[N]
            nndist_v[ijetB] = nndist_v[N]
            dijdist_v[ijetB] = dijdist_v[N]
            nxv[ijetB] = nxv[N]
            nyv[ijetB] = nyv[N]
            nzv[ijetB] = nzv[N]
            E2pv[ijetB] = E2pv[N]
            E2p_scaled[ijetB] = E2p_scaled[N]
            beam_term[ijetB] = beam_term[N]
        end
        N -= 1
        @inbounds for i in 1:N
            if (ijetB != N + 1) && (nni_v[i] == N + 1)
                nni_v[i] = ijetB
            else
                if (nni_v[i] == ijetA) || (nni_v[i] == ijetB) || (nni_v[i] > N)
                    update_nn_no_cross_arrays_precomputed!(nndist_v, nni_v, nxv, nyv, nzv,
                                                           E2p_scaled, beam_term, dijdist_v,
                                                           i, N, dij_factor, p, γ, R)
                end
            end
        end
        if ijetA != ijetB
            update_nn_cross_arrays_precomputed!(nndist_v, nni_v, nxv, nyv, nzv,
                                                E2p_scaled, beam_term, dijdist_v,
                                                ijetA, N, dij_factor, p, γ, R)
        end
    end
    clusterseq
end
