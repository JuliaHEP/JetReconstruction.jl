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
Base.@propagate_inbounds @inline function angular_distance(eereco, i, j)
    if hasproperty(eereco, :nx)
        nx = eereco.nx; ny = eereco.ny; nz = eereco.nz
        @muladd 1.0 - nx[i] * nx[j] - ny[i] * ny[j] - nz[i] * nz[j]
    else
        # Fallback for Array-of-structs (AoS)
        @muladd 1.0 - eereco[i].nx * eereco[j].nx - eereco[i].ny * eereco[j].ny - eereco[i].nz * eereco[j].nz
    end
end

"""
    valencia_distance(eereco, i, j, R) -> Float64

Calculate the Valencia distance between two jets `i` and `j` as
``min(E_i^{2β}, E_j^{2β}) * 2 * (1 - cos(θ_{ij})) / R²``.

# Arguments
- `eereco`: The array of `EERecoJet` objects.
- `i`: The first jet.
- `j`: The second jet.
- `R`: The jet radius parameter.

# Returns
- `Float64`: The Valencia distance between `i` and `j`.
"""
Base.@propagate_inbounds @inline function valencia_distance_inv(eereco, i, j, invR2)
    if hasproperty(eereco, :nx)
        nx = eereco.nx; ny = eereco.ny; nz = eereco.nz; E2p = eereco.E2p
        angular_dist = angular_distance_arrays(nx, ny, nz, i, j)
        # Valencia dij : min(E_i^{2β}, E_j^{2β}) * 2 * (1 - cos θ) * invR2
        min(E2p[i], E2p[j]) * 2 * angular_dist * invR2
    else
        # Fallback for Array-of-structs
        angular_dist = 1.0 - eereco[i].nx * eereco[j].nx - eereco[i].ny * eereco[j].ny - eereco[i].nz * eereco[j].nz
        min(eereco[i].E2p, eereco[j].E2p) * 2 * angular_dist * invR2
    end
end

Base.@propagate_inbounds @inline function valencia_distance(eereco, i, j, R)
    return valencia_distance_inv(eereco, i, j, inv(R * R))
end

# Array-based helpers: operate directly on field vectors from StructArray to avoid
# repeated eereco[i] indexing which can be slower.
Base.@propagate_inbounds @inline function angular_distance_arrays(nx, ny, nz, i, j)
    @muladd 1.0 - nx[i] * nx[j] - ny[i] * ny[j] - nz[i] * nz[j]
end

Base.@propagate_inbounds @inline function valencia_distance_inv_arrays(E2p, nx, ny, nz, i, j, invR2)
    angular_dist = angular_distance_arrays(nx, ny, nz, i, j)
    min(E2p[i], E2p[j]) * 2 * angular_dist * invR2
end

# Scaled variant: accepts pre-multiplied E2p_scaled = E2p * (2 * invR2)
Base.@propagate_inbounds @inline function valencia_distance_inv_scaled_arrays(E2p_scaled, nx, ny, nz, i, j)
    angular_dist = angular_distance_arrays(nx, ny, nz, i, j)
    min(E2p_scaled[i], E2p_scaled[j]) * angular_dist
end

"""
    valencia_beam_distance(eereco, i, γ, β) -> Float64

Calculate the Valencia beam distance for jet `i` using the FastJet ValenciaPlugin
definition: ``d_iB = E_i^{2β} * (sin θ_i)^{2γ}``, where ``cos θ_i = nz``
for unit direction cosines. Since ``sin^2 θ = 1 - nz^2``, we implement
``d_iB = E_i^{2β} * (1 - nz^2)^γ``.

# Arguments
- `eereco`: The array of `EERecoJet` objects.
- `i`: The jet index.
- `γ`: The angular exponent parameter used in the Valencia beam distance.
- `β`: The energy exponent (same as `p` in our implementation).

# Returns
- `Float64`: The Valencia beam distance for jet `i`.
"""
Base.@propagate_inbounds @inline function valencia_beam_distance(eereco, i, γ, β)
    if hasproperty(eereco, :nz)
        nzv = eereco.nz; E2pv = eereco.E2p
        nz_i = nzv[i]
        sin2 = 1 - nz_i * nz_i
        E2p = E2pv[i]
    else
        nz_i = eereco[i].nz
        sin2 = 1 - nz_i * nz_i
        E2p = eereco[i].E2p
    end
    # Fast-paths for common γ values to avoid pow in hot loop
    if γ == 1.0
        return E2p * sin2
    elseif γ == 2.0
        return E2p * (sin2 * sin2)
    else
        return E2p * sin2^γ
    end
end

# Array-based helper for Valencia beam distance to avoid StructArray getindex in hot loops
Base.@propagate_inbounds @inline function valencia_beam_distance_arrays(E2p, nz, i, γ, β)
    nz_i = nz[i]
    sin2 = 1 - nz_i * nz_i
    E2p_i = E2p[i]
    # Fast-paths for common γ values to avoid pow in hot loop
    if γ == 1.0
        return E2p_i * sin2
    elseif γ == 2.0
        return E2p_i * (sin2 * sin2)
    else
        return E2p_i * sin2^γ
    end
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
    j == 0 && return large_dij
    @inbounds min(eereco[i].E2p, eereco[j].E2p) * dij_factor * eereco[i].nndist
end

@inline function dij_dist(eereco, i, j, dij_factor, ::Val{JetAlgorithm.EEKt}, R = 4.0)
    j == 0 && return large_dij
    @inbounds min(eereco[i].E2p, eereco[j].E2p) * dij_factor * eereco[i].nndist
end

"""
    dij_dist(eereco, i, j, dij_factor, ::Val{JetAlgorithm.Valencia}, R)

Valencia dij distance uses the full Valencia metric, including the 2*(1-cosθ)/R² factor.
"""
@inline function dij_dist(eereco, i, j, dij_factor, ::Val{JetAlgorithm.Valencia}, R)
    j == 0 && return large_dij
    @inbounds valencia_distance(eereco, i, j, R)
end

# Fallback if a non-Algorithm token is passed
@inline function dij_dist(eereco, i, j, dij_factor, algorithm, R = 4.0)
    throw(ArgumentError("Algorithm $algorithm not supported for dij_dist"))
end

function get_angular_nearest_neighbours!(eereco, algorithm, dij_factor, p, γ = 1.0, R = 4.0)
    # Get the initial nearest neighbours for each jet
    N = length(eereco)
    # For Valencia, nearest-neighbour must be chosen on the full dij metric (FastJet NNH behaviour)
        @inbounds for i in 1:N
            local_nndist_i = Inf
            local_nni_i = i
            eereco.nndist[i] = local_nndist_i
            eereco.nni[i] = local_nni_i
        end
    # Nearest neighbour search
    @inbounds for i in 1:N
        @inbounds for j in (i + 1):N
            # Metric used to pick the nearest neighbour
            if algorithm == JetAlgorithm.Valencia
                # Use array helpers to avoid repeated StructArray getindex
                this_metric = valencia_distance_inv_arrays(eereco.E2p, eereco.nx, eereco.ny, eereco.nz, i, j, inv(R * R))
            else
                this_metric = angular_distance_arrays(eereco.nx, eereco.ny, eereco.nz, i, j)
            end

            # Using these ternary operators is faster than the if-else block
            better_nndist_i = this_metric < eereco[i].nndist
            eereco.nndist[i] = better_nndist_i ? this_metric : eereco.nndist[i]
            eereco.nni[i] = better_nndist_i ? j : eereco.nni[i]
            better_nndist_j = this_metric < eereco[j].nndist
            eereco.nndist[j] = better_nndist_j ? this_metric : eereco.nndist[j]
            eereco.nni[j] = better_nndist_j ? i : eereco.nni[j]
        end
    end
    # Nearest neighbour dij distance
    @inbounds for i in 1:N
        if algorithm == JetAlgorithm.Valencia
            eereco.dijdist[i] = valencia_distance(eereco, i, eereco[i].nni, R)
        else
            eereco.dijdist[i] = dij_dist(eereco, i, eereco[i].nni, dij_factor, algorithm, R)
        end
    end
    # For the EEKt algorithm, we need to check the beam distance as well
    # (This is structured to only check for EEKt once)
    if algorithm == JetAlgorithm.EEKt
        @inbounds for i in 1:N
            beam_closer = eereco[i].E2p < eereco[i].dijdist
            eereco.dijdist[i] = beam_closer ? eereco[i].E2p : eereco.dijdist[i]
            eereco.nni[i] = beam_closer ? 0 : eereco.nni[i]
        end
    elseif algorithm == JetAlgorithm.Valencia
        # Use array-based helper to avoid StructArray property checks and
        # reduce per-iteration overhead in the hot loop.
        E2p = eereco.E2p; nz = eereco.nz
        @inbounds for i in 1:N
            valencia_beam_dist = valencia_beam_distance_arrays(E2p, nz, i, γ, p)
            beam_closer = valencia_beam_dist < eereco.dijdist[i]
            eereco.dijdist[i] = beam_closer ? valencia_beam_dist : eereco.dijdist[i]
            eereco.nni[i] = beam_closer ? 0 : eereco.nni[i]
        end
    end
end

# Val-specialized nearest neighbour search (removes runtime branches in hot loops)
@inline function get_angular_nearest_neighbours!(eereco, ::Val{JetAlgorithm.Durham},
                                                 dij_factor, p, γ = 1.0, R = 4.0)
    N = length(eereco)
    nx = eereco.nx; ny = eereco.ny; nz = eereco.nz
    nndist = eereco.nndist; nni = eereco.nni
    @inbounds for i in 1:N
        local_nndist_i = large_distance
        local_nni_i = i
        @inbounds for j in 1:N
            if j != i
                this_metric = angular_distance_arrays(nx, ny, nz, i, j)
                better_nndist_i = this_metric < local_nndist_i
                local_nndist_i = better_nndist_i ? this_metric : local_nndist_i
                local_nni_i = better_nndist_i ? j : local_nni_i
                better_nndist_j = this_metric < nndist[j]
                nndist[j] = better_nndist_j ? this_metric : nndist[j]
                nni[j] = better_nndist_j ? i : nni[j]
            end
        end
        nndist[i] = local_nndist_i
        nni[i] = local_nni_i
    end
    @inbounds for i in 1:N
        eereco.dijdist[i] = dij_dist(eereco, i, nni[i], dij_factor,
                                     Val(JetAlgorithm.Durham), R)
    end
end

@inline function get_angular_nearest_neighbours!(eereco, ::Val{JetAlgorithm.EEKt},
                                                 dij_factor, p, γ = 1.0, R = 4.0)
    N = length(eereco)
    @inbounds for i in 1:N
        @inbounds for j in (i + 1):N
            this_metric = angular_distance(eereco, i, j)
            better_nndist_i = this_metric < eereco[i].nndist
            eereco.nndist[i] = better_nndist_i ? this_metric : eereco.nndist[i]
            eereco.nni[i] = better_nndist_i ? j : eereco.nni[i]
            better_nndist_j = this_metric < eereco[j].nndist
            eereco.nndist[j] = better_nndist_j ? this_metric : eereco.nndist[j]
            eereco.nni[j] = better_nndist_j ? i : eereco.nni[j]
        end
    end
    @inbounds for i in 1:N
        eereco.dijdist[i] = dij_dist(eereco, i, eereco[i].nni, dij_factor,
                                     Val(JetAlgorithm.EEKt), R)
    end
    @inbounds for i in 1:N
        beam_closer = eereco[i].E2p < eereco[i].dijdist
        eereco.dijdist[i] = beam_closer ? eereco[i].E2p : eereco.dijdist[i]
        eereco.nni[i] = beam_closer ? 0 : eereco.nni[i]
    end
end

@inline function get_angular_nearest_neighbours!(eereco, ::Val{JetAlgorithm.Valencia},
                                                 dij_factor, p, γ = 1.0, R = 4.0)
    # Fallback Val-specialised implementation kept for non-precomputed use.
    # The Valencia entrypoint uses a precomputed path (see _ee_genkt_algorithm_valencia)
    N = length(eereco)
    E2p = eereco.E2p; nx = eereco.nx; ny = eereco.ny; nz = eereco.nz
    nndist = eereco.nndist; nni = eereco.nni
    invR2 = inv(R * R)
        @inbounds for i in 1:N
            nndist[i] = Inf
            nni[i] = i
    end
    @inbounds for i in 1:N
        local_nndist_i = nndist[i]
        local_nni_i = nni[i]
        @inbounds for j in (i + 1):N
            this_metric = valencia_distance_inv_arrays(E2p, nx, ny, nz, i, j, invR2)
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
        eereco.dijdist[i] = valencia_distance_inv_arrays(E2p, nx, ny, nz, i, nni[i], invR2)
    end
    @inbounds for i in 1:N
        valencia_beam_dist = valencia_beam_distance(eereco, i, γ, p)
        beam_closer = valencia_beam_dist < eereco[i].dijdist
        eereco.dijdist[i] = beam_closer ? valencia_beam_dist : eereco.dijdist[i]
        eereco.nni[i] = beam_closer ? 0 : eereco.nni[i]
    end
end

## Precomputed Valencia nearest-neighbour initializer using precomputed arrays
Base.@propagate_inbounds @inline function get_angular_nearest_neighbours_valencia_precomputed!(eereco,
                                                                                              E2p_scaled::AbstractVector,
                                                                                              beam_term::AbstractVector,
                                                                                              p, γ = 1.0, R = 4.0)
    N = length(eereco)
    nx = eereco.nx; ny = eereco.ny; nz = eereco.nz
    nndist = eereco.nndist; nni = eereco.nni
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
        eereco.dijdist[i] = valencia_distance_inv_scaled_arrays(E2p_scaled, nx, ny, nz, i, nni[i])
    end
    @inbounds for i in 1:N
        beam_closer = beam_term[i] < eereco.dijdist[i]
        eereco.dijdist[i] = beam_closer ? beam_term[i] : eereco.dijdist[i]
        eereco.nni[i] = beam_closer ? 0 : eereco.nni[i]
    end
end

# Forwarder to Val-specialized version
@inline function get_angular_nearest_neighbours!(eereco,
                                                 algorithm::JetAlgorithm.Algorithm,
                                                 dij_factor, p, γ = 1.0, R = 4.0)
    return get_angular_nearest_neighbours!(eereco, Val(algorithm), dij_factor, p, γ, R)
end

# Update the nearest neighbour for jet i, w.r.t. all other active jets
@inline function update_nn_no_cross!(eereco, i, N, algorithm::JetAlgorithm.Algorithm,
                                     dij_factor, β = 1.0, γ = 1.0, R = 4.0)
    # Forward to Val-specialized implementations to avoid runtime branches
    return update_nn_no_cross!(eereco, i, N, Val(algorithm), dij_factor, β, γ, R)
end

# Val-specialized no-cross update
@inline function update_nn_no_cross!(eereco, i, N, ::Val{JetAlgorithm.Durham}, dij_factor,
                                     β = 1.0, γ = 1.0, R = 4.0)
    nndist = eereco.nndist; nni = eereco.nni
    nx = eereco.nx; ny = eereco.ny; nz = eereco.nz
    nndist[i] = large_distance
    nni[i] = i
    @inbounds for j in 1:(i-1)
        this_metric = angular_distance_arrays(nx, ny, nz, i, j)
        better_nndist_i = this_metric < nndist[i]
        nndist[i] = better_nndist_i ? this_metric : nndist[i]
        nni[i] = better_nndist_i ? j : nni[i]
    end
    @inbounds for j in (i+1):N
        this_metric = angular_distance_arrays(nx, ny, nz, i, j)
        better_nndist_i = this_metric < nndist[i]
        nndist[i] = better_nndist_i ? this_metric : nndist[i]
        nni[i] = better_nndist_i ? j : nni[i]
    end
    # Inline Durham dij computation to avoid function-call overhead
    E2p = eereco.E2p
    E2p_i = E2p[i]
    E2p_nni = E2p[nni[i]]
    minE2p = E2p_i < E2p_nni ? E2p_i : E2p_nni
    eereco.dijdist[i] = minE2p * dij_factor * nndist[i]
end

@inline function update_nn_no_cross!(eereco, i, N, ::Val{JetAlgorithm.EEKt}, dij_factor,
                                     β = 1.0, γ = 1.0, R = 4.0)
    nndist = eereco.nndist; nni = eereco.nni
    nx = eereco.nx; ny = eereco.ny; nz = eereco.nz
    E2p = eereco.E2p
    E2p_i = E2p[i]
    nndist[i] = large_distance
    nni[i] = i
    @inbounds for j in 1:(i-1)
        this_metric = angular_distance_arrays(nx, ny, nz, i, j)
        better_nndist_i = this_metric < nndist[i]
        nndist[i] = better_nndist_i ? this_metric : nndist[i]
        nni[i] = better_nndist_i ? j : nni[i]
    end
    @inbounds for j in (i+1):N
        this_metric = angular_distance_arrays(nx, ny, nz, i, j)
        better_nndist_i = this_metric < nndist[i]
        nndist[i] = better_nndist_i ? this_metric : nndist[i]
        nni[i] = better_nndist_i ? j : nni[i]
    end
    # Inline EEKt dij computation and beam check
    E2p_nni = E2p[nni[i]]
    minE2p = E2p_i < E2p_nni ? E2p_i : E2p_nni
    eereco.dijdist[i] = minE2p * dij_factor * nndist[i]
    beam_close = E2p_i < eereco.dijdist[i]
    eereco.dijdist[i] = beam_close ? E2p_i : eereco.dijdist[i]
    eereco.nni[i] = beam_close ? 0 : eereco.nni[i]
end

@inline function update_nn_no_cross!(eereco, i, N, ::Val{JetAlgorithm.Valencia},
                                     dij_factor, β = 1.0, γ = 1.0, R = 4.0)
    nndist = eereco.nndist; nni = eereco.nni
    nx = eereco.nx; ny = eereco.ny; nz = eereco.nz
    E2p = eereco.E2p
    E2p_i = E2p[i]
    nndist[i] = Inf
    nni[i] = i
    invR2 = inv(R * R)
    @inbounds for j in 1:(i-1)
        this_metric = valencia_distance_inv_arrays(E2p, nx, ny, nz, i, j, invR2)
        better_nndist_i = this_metric < nndist[i]
        nndist[i] = better_nndist_i ? this_metric : nndist[i]
        nni[i] = better_nndist_i ? j : nni[i]
    end
    @inbounds for j in (i+1):N
        this_metric = valencia_distance_inv_arrays(E2p, nx, ny, nz, i, j, invR2)
        better_nndist_i = this_metric < nndist[i]
        nndist[i] = better_nndist_i ? this_metric : nndist[i]
        nni[i] = better_nndist_i ? j : nni[i]
    end
    eereco.dijdist[i] = valencia_distance_inv_arrays(E2p, nx, ny, nz, i, nni[i], invR2)
    valencia_beam_dist = valencia_beam_distance(eereco, i, γ, β)
    beam_close = valencia_beam_dist < eereco.dijdist[i]
    eereco.dijdist[i] = beam_close ? valencia_beam_dist : eereco.dijdist[i]
    eereco.nni[i] = beam_close ? 0 : eereco.nni[i]
end

# (Forwarding handled earlier; avoid duplicate definition)

@inline function update_nn_cross!(eereco, i, N, algorithm::JetAlgorithm.Algorithm,
                                  dij_factor, β = 1.0, γ = 1.0, R = 4.0)
    # Forward to Val-specialized implementations to avoid runtime branches
    return update_nn_cross!(eereco, i, N, Val(algorithm), dij_factor, β, γ, R)
end

# Val-specialized cross update
@inline function update_nn_cross!(eereco, i, N, ::Val{JetAlgorithm.Durham}, dij_factor,
                                  β = 1.0, γ = 1.0, R = 4.0)
    nndist = eereco.nndist; nni = eereco.nni
    nx = eereco.nx; ny = eereco.ny; nz = eereco.nz
    nndist[i] = large_distance
    nni[i] = i
    E2p = eereco.E2p
    E2p_i = E2p[i]
    @inbounds for j in 1:N
        if j != i
            this_metric = angular_distance_arrays(nx, ny, nz, i, j)
            better_nndist_i = this_metric < nndist[i]
            nndist[i] = better_nndist_i ? this_metric : nndist[i]
            nni[i] = better_nndist_i ? j : nni[i]
            if this_metric < nndist[j]
                nndist[j] = this_metric
                nni[j] = i
                # Hoist E2p_j and compute new dij once
                E2p_j = E2p[j]
                new_dij = (E2p_i < E2p_j ? E2p_i : E2p_j) * dij_factor * this_metric
                eereco.dijdist[j] = new_dij
            end
        end
    end
    # Inline dij for i
    E2p_nni = E2p[nni[i]]
    minE2p_i = E2p_i < E2p_nni ? E2p_i : E2p_nni
    eereco.dijdist[i] = minE2p_i * dij_factor * nndist[i]
end

@inline function update_nn_cross!(eereco, i, N, ::Val{JetAlgorithm.EEKt}, dij_factor,
                                  β = 1.0, γ = 1.0, R = 4.0)
    nndist = eereco.nndist; nni = eereco.nni
    nx = eereco.nx; ny = eereco.ny; nz = eereco.nz
    E2p = eereco.E2p
    E2p_i = E2p[i]
    nndist[i] = large_distance
    nni[i] = i
    @inbounds for j in 1:N
        if j != i
            this_metric = angular_distance_arrays(nx, ny, nz, i, j)
            better_nndist_i = this_metric < nndist[i]
            nndist[i] = better_nndist_i ? this_metric : nndist[i]
            nni[i] = better_nndist_i ? j : nni[i]
            if this_metric < nndist[j]
                nndist[j] = this_metric
                nni[j] = i
                # Hoist E2p_j and compute new dij once, then beam-check
                E2p_j = E2p[j]
                new_dij = (E2p_i < E2p_j ? E2p_i : E2p_j) * dij_factor * this_metric
                if E2p_j < new_dij
                    eereco.dijdist[j] = E2p_j
                    nni[j] = 0
                else
                    eereco.dijdist[j] = new_dij
                end
            end
        end
    end
    # Inline dij for i and beam check
    E2p_nni = E2p[nni[i]]
    minE2p_i = E2p_i < E2p_nni ? E2p_i : E2p_nni
    eereco.dijdist[i] = minE2p_i * dij_factor * nndist[i]
    beam_close = E2p_i < eereco.dijdist[i]
    eereco.dijdist[i] = beam_close ? E2p_i : eereco.dijdist[i]
    nni[i] = beam_close ? 0 : nni[i]
end

@inline function update_nn_cross!(eereco, i, N, ::Val{JetAlgorithm.Valencia}, dij_factor,
                                  β = 1.0, γ = 1.0, R = 4.0)
    nndist = eereco.nndist; nni = eereco.nni
    nx = eereco.nx; ny = eereco.ny; nz = eereco.nz
    E2p = eereco.E2p
    nndist[i] = Inf
    nni[i] = i
    invR2 = inv(R * R)
    @inbounds for j in 1:N
        if j != i
            this_metric = valencia_distance_inv_arrays(E2p, nx, ny, nz, i, j, invR2)
            better_nndist_i = this_metric < nndist[i]
            nndist[i] = better_nndist_i ? this_metric : nndist[i]
            nni[i] = better_nndist_i ? j : nni[i]
            if this_metric < nndist[j]
                nndist[j] = this_metric
                nni[j] = i
                # Use the already-computed metric for dij (Valencia uses full dij here)
                eereco.dijdist[j] = this_metric
                valencia_beam_dist = valencia_beam_distance_arrays(E2p, nz, j, γ, β)
                if valencia_beam_dist < eereco.dijdist[j]
                    eereco.dijdist[j] = valencia_beam_dist
                    nni[j] = 0
                end
            end
        end
    end
    eereco.dijdist[i] = valencia_distance_inv_arrays(E2p, nx, ny, nz, i, nni[i], invR2)
    valencia_beam_dist = valencia_beam_distance_arrays(E2p, nz, i, γ, β)
    beam_close = valencia_beam_dist < eereco.dijdist[i]
    eereco.dijdist[i] = beam_close ? valencia_beam_dist : eereco.dijdist[i]
    nni[i] = beam_close ? 0 : nni[i]
end

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

################################################################################
# Array-based nearest-neighbour update helpers (operate on raw vectors)
################################################################################

@inline function update_nn_no_cross_arrays!(nndist::AbstractVector, nni::AbstractVector,
                                           nx::AbstractVector, ny::AbstractVector, nz::AbstractVector,
                                           E2p::AbstractVector, dijdist::AbstractVector,
                                           i::Integer, N::Integer, algorithm::JetAlgorithm.Algorithm,
                                           dij_factor, β = 1.0, γ = 1.0, R = 4.0)
    return update_nn_no_cross_arrays!(nndist, nni, nx, ny, nz, E2p, dijdist, i, N, Val(algorithm), dij_factor, β, γ, R)
end

@inline function update_nn_no_cross_arrays!(nndist::AbstractVector, nni::AbstractVector,
                                           nx::AbstractVector, ny::AbstractVector, nz::AbstractVector,
                                           E2p::AbstractVector, dijdist::AbstractVector,
                                           i::Integer, N::Integer, ::Val{JetAlgorithm.Durham},
                                           dij_factor, β = 1.0, γ = 1.0, R = 4.0)
    nndist[i] = large_distance
    nni[i] = i
    @inbounds for j in 1:N
        if j != i
            this_metric = angular_distance_arrays(nx, ny, nz, i, j)
            better_nndist_i = this_metric < nndist[i]
            nndist[i] = better_nndist_i ? this_metric : nndist[i]
            nni[i] = better_nndist_i ? j : nni[i]
        end
    end
    # Inline Durham dij computation
    E2p_i = E2p[i]
    E2p_nni = E2p[nni[i]]
    minE2p = E2p_i < E2p_nni ? E2p_i : E2p_nni
    dijdist[i] = minE2p * dij_factor * nndist[i]
end

@inline function update_nn_no_cross_arrays!(nndist::AbstractVector, nni::AbstractVector,
                                           nx::AbstractVector, ny::AbstractVector, nz::AbstractVector,
                                           E2p::AbstractVector, dijdist::AbstractVector,
                                           i::Integer, N::Integer, ::Val{JetAlgorithm.EEKt},
                                           dij_factor, β = 1.0, γ = 1.0, R = 4.0)
    nndist[i] = large_distance
    nni[i] = i
    E2p_i = E2p[i]
    @inbounds for j in 1:N
        if j != i
            this_metric = angular_distance_arrays(nx, ny, nz, i, j)
            better_nndist_i = this_metric < nndist[i]
            nndist[i] = better_nndist_i ? this_metric : nndist[i]
            nni[i] = better_nndist_i ? j : nni[i]
        end
    end
    E2p_nni = E2p[nni[i]]
    minE2p = E2p_i < E2p_nni ? E2p_i : E2p_nni
    dijdist[i] = minE2p * dij_factor * nndist[i]
    beam_close = E2p_i < dijdist[i]
    dijdist[i] = beam_close ? E2p_i : dijdist[i]
    nni[i] = beam_close ? 0 : nni[i]
end

@inline function update_nn_no_cross_arrays!(nndist::AbstractVector, nni::AbstractVector,
                                           nx::AbstractVector, ny::AbstractVector, nz::AbstractVector,
                                           E2p::AbstractVector, dijdist::AbstractVector,
                                           i::Integer, N::Integer, ::Val{JetAlgorithm.Valencia},
                                           dij_factor, β = 1.0, γ = 1.0, R = 4.0)
    # Precomputed variant lives in a separate function; keep this fallback here.
    invR2 = inv(R * R)
    nndist_i = Inf
    nni_i = i
    @inbounds for j in 1:N
        if j != i
            this_metric = valencia_distance_inv_arrays(E2p, nx, ny, nz, i, j, invR2)
            if this_metric < nndist_i
                nndist_i = this_metric
                nni_i = j
            end
        end
    end
    # compute dijdist and beam check using locals then write back once
    dijdist_i = valencia_distance_inv_arrays(E2p, nx, ny, nz, i, nni_i, invR2)
    valencia_beam_dist = valencia_beam_distance_arrays(E2p, nz, i, γ, β)
    if valencia_beam_dist < dijdist_i
        dijdist_i = valencia_beam_dist
        nni_i = 0
    end
    nndist[i] = nndist_i
    nni[i] = nni_i
    dijdist[i] = dijdist_i
end

## Precomputed Valencia no-cross update using E2p_scaled and beam_term
@inline function update_nn_no_cross_arrays_precomputed!(nndist::AbstractVector, nni::AbstractVector,
                                                       nx::AbstractVector, ny::AbstractVector, nz::AbstractVector,
                                                       E2p_scaled::AbstractVector, beam_term::AbstractVector,
                                                       dijdist::AbstractVector,
                                                       i::Integer, N::Integer,
                                                       dij_factor, β = 1.0, γ = 1.0, R = 4.0)
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

@inline function update_nn_cross_arrays!(nndist::AbstractVector, nni::AbstractVector,
                                        nx::AbstractVector, ny::AbstractVector, nz::AbstractVector,
                                        E2p::AbstractVector, dijdist::AbstractVector,
                                        i::Integer, N::Integer, algorithm::JetAlgorithm.Algorithm,
                                        dij_factor, β = 1.0, γ = 1.0, R = 4.0)
    return update_nn_cross_arrays!(nndist, nni, nx, ny, nz, E2p, dijdist, i, N, Val(algorithm), dij_factor, β, γ, R)
end

@inline function update_nn_cross_arrays!(nndist::AbstractVector, nni::AbstractVector,
                                        nx::AbstractVector, ny::AbstractVector, nz::AbstractVector,
                                        E2p::AbstractVector, dijdist::AbstractVector,
                                        i::Integer, N::Integer, ::Val{JetAlgorithm.Durham},
                                        dij_factor, β = 1.0, γ = 1.0, R = 4.0)
    nndist[i] = large_distance
    nni[i] = i
    E2p_i = E2p[i]
    @inbounds for j in 1:(i-1)
        this_metric = angular_distance_arrays(nx, ny, nz, i, j)
        better_nndist_i = this_metric < nndist[i]
        nndist[i] = better_nndist_i ? this_metric : nndist[i]
        nni[i] = better_nndist_i ? j : nni[i]
        if this_metric < nndist[j]
            nndist[j] = this_metric
            nni[j] = i
            E2p_j = E2p[j]
            dijdist[j] = (E2p_i < E2p_j ? E2p_i : E2p_j) * dij_factor * this_metric
        end
    end
    @inbounds for j in (i+1):N
        this_metric = angular_distance_arrays(nx, ny, nz, i, j)
        better_nndist_i = this_metric < nndist[i]
        nndist[i] = better_nndist_i ? this_metric : nndist[i]
        nni[i] = better_nndist_i ? j : nni[i]
        if this_metric < nndist[j]
            nndist[j] = this_metric
            nni[j] = i
            E2p_j = E2p[j]
            minE2p = E2p_i < E2p_j ? E2p_i : E2p_j
            dijdist[j] = minE2p * dij_factor * this_metric
        end
    end
    E2p_nni = E2p[nni[i]]
    minE2p_i = E2p_i < E2p_nni ? E2p_i : E2p_nni
    dijdist[i] = minE2p_i * dij_factor * nndist[i]
end

@inline function update_nn_cross_arrays!(nndist::AbstractVector, nni::AbstractVector,
                                        nx::AbstractVector, ny::AbstractVector, nz::AbstractVector,
                                        E2p::AbstractVector, dijdist::AbstractVector,
                                        i::Integer, N::Integer, ::Val{JetAlgorithm.EEKt},
                                        dij_factor, β = 1.0, γ = 1.0, R = 4.0)
    E2p_i = E2p[i]
    nndist[i] = large_distance
    nni[i] = i
    @inbounds for j in 1:(i-1)
        this_metric = angular_distance_arrays(nx, ny, nz, i, j)
        better_nndist_i = this_metric < nndist[i]
        nndist[i] = better_nndist_i ? this_metric : nndist[i]
        nni[i] = better_nndist_i ? j : nni[i]
            if this_metric < nndist[j]
                nndist[j] = this_metric
                nni[j] = i
                E2p_j = E2p[j]
                new_dij = (E2p_i < E2p_j ? E2p_i : E2p_j) * dij_factor * this_metric
                if E2p_j < new_dij
                    dijdist[j] = E2p_j
                    nni[j] = 0
                else
                    dijdist[j] = new_dij
                end
            end
    end
    @inbounds for j in (i+1):N
        this_metric = angular_distance_arrays(nx, ny, nz, i, j)
        better_nndist_i = this_metric < nndist[i]
        nndist[i] = better_nndist_i ? this_metric : nndist[i]
        nni[i] = better_nndist_i ? j : nni[i]
        if this_metric < nndist[j]
            nndist[j] = this_metric
            nni[j] = i
            E2p_j = E2p[j]
            minE2p = E2p_i < E2p_j ? E2p_i : E2p_j
            dijdist[j] = minE2p * dij_factor * this_metric
            if E2p_j < dijdist[j]
                dijdist[j] = E2p_j
                nni[j] = 0
            end
        end
    end
    E2p_nni = E2p[nni[i]]
    minE2p_i = E2p_i < E2p_nni ? E2p_i : E2p_nni
    dijdist[i] = minE2p_i * dij_factor * nndist[i]
    beam_close = E2p_i < dijdist[i]
    dijdist[i] = beam_close ? E2p_i : dijdist[i]
    nni[i] = beam_close ? 0 : nni[i]
end

@inline function update_nn_cross_arrays!(nndist::AbstractVector, nni::AbstractVector,
                                        nx::AbstractVector, ny::AbstractVector, nz::AbstractVector,
                                        E2p::AbstractVector, dijdist::AbstractVector,
                                        i::Integer, N::Integer, ::Val{JetAlgorithm.Valencia},
                                        dij_factor, β = 1.0, γ = 1.0, R = 4.0)
    # Operate on locals for slot i to reduce setindex traffic; updates to other
    # slots (j) still write directly since they modify different indices.
    # Precomputed variant fallback uses the non-precomputed helpers; there is
    # a separate precomputed cross-update below.
    invR2 = inv(R * R)
    nndist_i = Inf
    nni_i = i
    @inbounds for j in 1:(i-1)
        this_metric = valencia_distance_inv_arrays(E2p, nx, ny, nz, i, j, invR2)
        if this_metric < nndist_i
            nndist_i = this_metric
            nni_i = j
        end
        if this_metric < nndist[j]
            nndist[j] = this_metric
            nni[j] = i
            dijdist[j] = this_metric
            valencia_beam_dist = valencia_beam_distance_arrays(E2p, nz, j, γ, β)
            if valencia_beam_dist < dijdist[j]
                dijdist[j] = valencia_beam_dist
                nni[j] = 0
            end
        end
    end
    @inbounds for j in (i+1):N
        this_metric = valencia_distance_inv_arrays(E2p, nx, ny, nz, i, j, invR2)
        if this_metric < nndist_i
            nndist_i = this_metric
            nni_i = j
        end
        if this_metric < nndist[j]
            nndist[j] = this_metric
            nni[j] = i
            dijdist[j] = this_metric
            valencia_beam_dist = valencia_beam_distance_arrays(E2p, nz, j, γ, β)
            if valencia_beam_dist < dijdist[j]
                dijdist[j] = valencia_beam_dist
                nni[j] = 0
            end
        end
    end
    # Finalize slot i
    dijdist_i = valencia_distance_inv_arrays(E2p, nx, ny, nz, i, nni_i, invR2)
    valencia_beam_dist = valencia_beam_distance_arrays(E2p, nz, i, γ, β)
    if valencia_beam_dist < dijdist_i
        dijdist_i = valencia_beam_dist
        nni_i = 0
    end
    nndist[i] = nndist_i
    nni[i] = nni_i
    dijdist[i] = dijdist_i
end

## Precomputed Valencia cross-update using E2p_scaled and beam_term
@inline function update_nn_cross_arrays_precomputed!(nndist::AbstractVector, nni::AbstractVector,
                                                    nx::AbstractVector, ny::AbstractVector, nz::AbstractVector,
                                                    E2p_scaled::AbstractVector, beam_term::AbstractVector,
                                                    dijdist::AbstractVector,
                                                    i::Integer, N::Integer,
                                                    dij_factor, β = 1.0, γ = 1.0, R = 4.0)
    nndist_i = Inf
    nni_i = i
    @inbounds for j in 1:(i-1)
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
    @inbounds for j in (i+1):N
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

    # For Valencia, if β is provided, overwrite p
    if algorithm == JetAlgorithm.Valencia && β !== nothing
        p = β
    end

    # Check for consistency algorithm power
    p = get_algorithm_power(p = p, algorithm = algorithm)

    # Integer p if possible, i.e. if not running Valencia
    if algorithm != JetAlgorithm.Valencia
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

    # Now call the actual reconstruction method, tuned for our internal EDM
    # Dispatch once on the algorithm to call a small number of specialised
    # internal implementations. This avoids per-iteration branching/Val() in
    # the hot inner loops.
    if algorithm == JetAlgorithm.Valencia
        return _ee_genkt_algorithm_valencia(particles = recombination_particles, p = p, R = R,
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
# Durham-specialised implementation (optimized inner loops using array helpers)
################################################################################
function _ee_genkt_algorithm_durham(; particles::AbstractVector{EEJet},
                                    algorithm::JetAlgorithm.Algorithm, p::Real, R = 4.0,
                                    recombine = addjets, γ::Real = 1.0,
                                    beta::Union{Real, Nothing} = nothing)
    # Bounds
    N::Int = length(particles)

    R2 = R^2
    if algorithm == JetAlgorithm.Valencia && beta !== nothing
        p = beta
    end

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
    index = eereco.index; nni = eereco.nni; nndist = eereco.nndist
    dijdist = eereco.dijdist; nx = eereco.nx; ny = eereco.ny; nz = eereco.nz
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
                    update_nn_no_cross_arrays!(nndist, nni, nx, ny, nz, E2p, dijdist,
                                               i, N, Val(JetAlgorithm.Durham), dij_factor, p, γ, R)
                end
            end
        end

        if ijetA != ijetB
            update_nn_cross_arrays!(nndist, nni, nx, ny, nz, E2p, dijdist,
                                    ijetA, N, Val(JetAlgorithm.Durham), dij_factor, p, γ, R)
        end
    end

    clusterseq
end

################################################################################
# Valencia-specialised implementation
################################################################################
function _ee_genkt_algorithm_valencia(; particles::AbstractVector{EEJet},
                                      algorithm::JetAlgorithm.Algorithm, p::Real, R = 4.0,
                                      recombine = addjets, γ::Real = 1.0,
                                      beta::Union{Real, Nothing} = nothing)
    # Bounds
    N::Int = length(particles)

    R2 = R^2
    # Valencia uses p as β when passed through
    if algorithm == JetAlgorithm.Valencia && beta !== nothing
        p = beta
    end

    # Constant factor for the Valencia dij metric
    dij_factor = 1.0

    # For optimised reconstruction generate an SoA containing the necessary
    # jet information and populate it accordingly
    eereco = StructArray{EERecoJet}(undef, N)

    fill_reco_array!(eereco, particles, R2, p)

    # Setup the initial history and get the total energy
    history, Qtot = initial_history(particles)

    clusterseq = ClusterSequence(algorithm, p, R, RecoStrategy.N2Plain, particles, history,
                                 Qtot)

    # Run over initial pairs of jets to find nearest neighbours (precomputed path)
    # prepare precomputed arrays below before calling the initialized helper
    

    # Alias StructArray fields into local vectors to avoid allocations from
    # copying while still avoiding repeated StructArray field lookups in
    # hot loops. These locals point directly at the underlying vectors, so
    # no explicit writeback is required at the end.
    indexv = eereco.index
    nni_v = eereco.nni
    nndist_v = eereco.nndist
    dijdist_v = eereco.dijdist
    nxv = eereco.nx
    nyv = eereco.ny
    nzv = eereco.nz
    E2pv = eereco.E2p

    # Precompute scaled E2p and beam_term for Valencia to avoid repeated work
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

    # Now run NN init using precomputed helpers
    get_angular_nearest_neighbours_valencia_precomputed!(eereco, E2p_scaled, beam_term, p, γ, R)

    # Now we can start the main loop
    iter = 0
    while N != 0
        iter += 1

    dij_min, ijetA = fast_findmin(dijdist_v, N)
    ijetB = nni_v[ijetA]

        # Now we check if there is a "beam" merge possibility
            if ijetB == 0
            # Shouldn't happen for Valencia (beam handled via valencia_beam checks)
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

            # Insert merged jet into our local SoA (avoid writing StructArray)
            indexv[ijetA] = newjet_k
            nni_v[ijetA] = 0
            nndist_v[ijetA] = R2
            nxv[ijetA] = nx(merged_jet)
            nyv[ijetA] = ny(merged_jet)
            nzv[ijetA] = nz(merged_jet)
            # Compute E2p like insert_new_jet! would
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
            # Recompute precomputed derived arrays for the merged slot
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
            # Local copy from slot N -> slot ijetB (avoid StructArray ops)
            indexv[ijetB] = indexv[N]
            nni_v[ijetB] = nni_v[N]
            nndist_v[ijetB] = nndist_v[N]
            dijdist_v[ijetB] = dijdist_v[N]
            nxv[ijetB] = nxv[N]
            nyv[ijetB] = nyv[N]
            nzv[ijetB] = nzv[N]
            E2pv[ijetB] = E2pv[N]
            # Also copy precomputed derived arrays
            E2p_scaled[ijetB] = E2p_scaled[N]
            beam_term[ijetB] = beam_term[N]
        end

        N -= 1

    # Update nearest neighbours step
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

    # Locals alias the StructArray fields directly, so there is no separate
    # writeback step required here.

    clusterseq
end

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

    # If Valencia requested, forward to the Valencia-specialised path
    if algorithm == JetAlgorithm.Valencia
        return _ee_genkt_algorithm_valencia(particles = particles, algorithm = algorithm,
                                            p = p, R = R, recombine = recombine,
                                            γ = γ, beta = beta)
    end

    R2 = R^2
    if algorithm == JetAlgorithm.Valencia && beta !== nothing
        p = beta
    end

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
