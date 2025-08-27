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
    dij_dist(eereco, i, j, dij_factor, algorithm::JetAlgorithm.Algorithm, invR2)

Compute dij distance for (Durham, EEKt, Valencia) using simple conditionals.
Beam index `j==0` returns a large sentinel. For Valencia we use the full
Valencia metric (independent of `dij_factor`).
"""
@inline function dij_dist(eereco, i, j, dij_factor, algorithm, invR2)
    j == 0 && return large_dij
    @inbounds begin
        if algorithm == JetAlgorithm.Valencia
            return valencia_distance(eereco, i, j, invR2)
        else
            # Durham & EEKt share same form here (min(E2p_i,E2p_j) * dij_factor * angular_metric)
            return min(eereco[i].E2p, eereco[j].E2p) * dij_factor * eereco[i].nndist
        end
    end
end

"""
    valencia_distance(eereco, i, j, invR2) -> Float64

Calculate the Valencia distance between two jets `i` and `j` as
``min(E_i^{2β}, E_j^{2β}) * 2 * (1 - cos(θ_{ij})) * invR2``.

# Arguments
- `eereco`: The array of `EERecoJet` objects.
- `i`: The first jet.
- `j`: The second jet.
- `invR2`: The inverse square of the radius, i.e. ``1 / R^2``.

# Returns
- `Float64`: The Valencia distance between `i` and `j`.
"""
Base.@propagate_inbounds @inline function valencia_distance(eereco, i, j, invR2)
    @muladd angular_dist = 1.0 - eereco[i].nx * eereco[j].nx - eereco[i].ny * eereco[j].ny -
                           eereco[i].nz * eereco[j].nz
    return min(eereco[i].E2p, eereco[j].E2p) * 2 * angular_dist * invR2
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
    nz_i = eereco[i].nz
    sin2 = 1 - nz_i * nz_i
    E2p = eereco[i].E2p
    if γ == 1.0
        return E2p * sin2
    elseif γ == 2.0
        return E2p * (sin2 * sin2)
    else
        return E2p * sin2^γ
    end
end

function get_angular_nearest_neighbours!(eereco, algorithm, dij_factor, p, invR2, γ = 1.0)
    # Get the initial nearest neighbours for each jet
    N = length(eereco)
    # Initialise sentinels so the first comparison always wins
    @inbounds for i in 1:N
        if algorithm == JetAlgorithm.Valencia
            eereco.nndist[i] = Inf
        else
            eereco.nndist[i] = large_distance
        end
        eereco.nni[i] = i
    end
    # Nearest neighbour search
    @inbounds for i in 1:N
        @inbounds for j in (i + 1):N
            # Metric used to pick the nearest neighbour
            if algorithm == JetAlgorithm.Valencia
                # Use canonical Valencia distance (StructArray-aware)
                this_metric = valencia_distance(eereco, i, j, invR2)
            else
                this_metric = angular_distance(eereco, i, j)
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
            eereco.dijdist[i] = valencia_distance(eereco, i, eereco[i].nni, invR2)
        else
            eereco.dijdist[i] = dij_dist(eereco, i, eereco[i].nni, dij_factor, algorithm,
                                         invR2)
        end
    end
    # For the EEKt and Valencia algorithms, we need to check the beam distance as well
    # (This is structured to check each algorithm's beam distance once)
    if algorithm == JetAlgorithm.EEKt
        @inbounds for i in 1:N
            beam_closer = eereco[i].E2p < eereco[i].dijdist
            eereco.dijdist[i] = beam_closer ? eereco[i].E2p : eereco.dijdist[i]
            eereco.nni[i] = beam_closer ? 0 : eereco.nni[i]
        end
    elseif algorithm == JetAlgorithm.Valencia
        @inbounds for i in 1:N
            valencia_beam_dist = valencia_beam_distance(eereco, i, γ, p)
            beam_closer = valencia_beam_dist < eereco[i].dijdist
            eereco.dijdist[i] = beam_closer ? valencia_beam_dist : eereco.dijdist[i]
            eereco.nni[i] = beam_closer ? 0 : eereco.nni[i]
        end
    end
end

@inline function update_nn_no_cross!(eereco, i, N, algorithm::JetAlgorithm.Algorithm,
                                     dij_factor, invR2, β = 1.0, γ = 1.0)
    # Valencia metric is unbounded, others use a large finite value
    if algorithm == JetAlgorithm.Valencia
        eereco.nndist[i] = Inf
    else
        eereco.nndist[i] = large_distance
    end
    eereco.nni[i] = i
    # Scan all other jets, update i, and cross-update j
    @inbounds for j in 1:N
        if j != i
            this_metric = algorithm == JetAlgorithm.Valencia ?
                          valencia_distance(eereco, i, j, invR2) :
                          angular_distance(eereco, i, j)
            better_i = this_metric < eereco[i].nndist
            eereco.nndist[i] = better_i ? this_metric : eereco.nndist[i]
            eereco.nni[i] = better_i ? j : eereco.nni[i]
        end
    end
    # Set dij for i using unified dispatcher and apply beam checks
    eereco.dijdist[i] = dij_dist(eereco, i, eereco[i].nni, dij_factor, algorithm, invR2)
    if algorithm == JetAlgorithm.EEKt
        beam_close = eereco[i].E2p < eereco[i].dijdist
        eereco.dijdist[i] = beam_close ? eereco[i].E2p : eereco.dijdist[i]
        eereco.nni[i] = beam_close ? 0 : eereco.nni[i]
    elseif algorithm == JetAlgorithm.Valencia
        val_beam = valencia_beam_distance(eereco, i, γ, β)
        beam_close = val_beam < eereco[i].dijdist
        eereco.dijdist[i] = beam_close ? val_beam : eereco.dijdist[i]
        eereco.nni[i] = beam_close ? 0 : eereco.nni[i]
    end
end

@inline function update_nn_cross!(eereco, i, N, algorithm::JetAlgorithm.Algorithm,
                                  dij_factor, invR2, β = 1.0, γ = 1.0)
    # Valencia metric is unbounded, others use a large finite value
    if algorithm == JetAlgorithm.Valencia
        eereco.nndist[i] = Inf
    else
        eereco.nndist[i] = large_distance
    end
    eereco.nni[i] = i
    # Scan all other jets, update i, and cross-update j
    @inbounds for j in 1:N
        if j != i
            this_metric = algorithm == JetAlgorithm.Valencia ?
                          valencia_distance(eereco, i, j, invR2) :
                          angular_distance(eereco, i, j)
            better_i = this_metric < eereco[i].nndist
            eereco.nndist[i] = better_i ? this_metric : eereco.nndist[i]
            eereco.nni[i] = better_i ? j : eereco.nni[i]

            if this_metric < eereco[j].nndist
                eereco.nndist[j] = this_metric
                eereco.nni[j] = i
                eereco.dijdist[j] = dij_dist(eereco, j, i, dij_factor, algorithm, invR2)
                if algorithm == JetAlgorithm.EEKt
                    if eereco[j].E2p < eereco[j].dijdist
                        eereco.dijdist[j] = eereco[j].E2p
                        eereco.nni[j] = 0
                    end
                elseif algorithm == JetAlgorithm.Valencia
                    val_beam_j = valencia_beam_distance(eereco, j, γ, β)
                    if val_beam_j < eereco[j].dijdist
                        eereco.dijdist[j] = val_beam_j
                        eereco.nni[j] = 0
                    end
                end
            end
        end
    end
    # Set dij for i using unified dispatcher and apply beam checks
    eereco.dijdist[i] = dij_dist(eereco, i, eereco[i].nni, dij_factor, algorithm, invR2)
    if algorithm == JetAlgorithm.EEKt
        beam_close = eereco[i].E2p < eereco[i].dijdist
        eereco.dijdist[i] = beam_close ? eereco[i].E2p : eereco.dijdist[i]
        eereco.nni[i] = beam_close ? 0 : eereco.nni[i]
    elseif algorithm == JetAlgorithm.Valencia
        val_beam_i = valencia_beam_distance(eereco, i, γ, β)
        beam_close = val_beam_i < eereco[i].dijdist
        eereco.dijdist[i] = beam_close ? val_beam_i : eereco.dijdist[i]
        eereco.nni[i] = beam_close ? 0 : eereco.nni[i]
    end
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

Base.@propagate_inbounds @inline function fill_reco_array!(eereco, particles, invR2, p)
    @inbounds for i in eachindex(particles)
        eereco.index[i] = i
        eereco.nni[i] = 0
        eereco.nndist[i] = inv(invR2) # R^2 as initial sentinel for angular algorithms
        # eereco.dijdist[i] = UNDEF # Does not need to be initialised
        eereco.nx[i] = nx(particles[i])
        eereco.ny[i] = ny(particles[i])
        eereco.nz[i] = nz(particles[i])
        eereco.E2p[i] = energy(particles[i])^(2p)
        # No precomputed beam factor; compute on demand to preserve previous behaviour
    end
end

Base.@propagate_inbounds @inline function insert_new_jet!(eereco, i, newjet_k, invR2,
                                                          merged_jet, p)
    @inbounds begin
        eereco.index[i] = newjet_k
        eereco.nni[i] = 0
        eereco.nndist[i] = inv(invR2)
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
                       preprocess = nothing, γ::Real = 1.0, β::Union{Real, Nothing} = nothing) where {T}

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
- `β::Union{Real, Nothing} = nothing`: Optional alias for the Valencia energy exponent; if provided for
    Valencia it overrides `p`.

# Returns
- The result of the jet clustering as a `ClusterSequence` object.

# Notes
This is the public interface to the e+e- jet clustering algorithm. The function
will check for consistency between the algorithm and the power parameter as
needed. It will then prepare the internal EDM particles for the clustering
itself, and call the actual reconstruction method `_ee_genkt_algorithm!`.

If the algorithm is Durham, `R` is nominally set to 4.
If the algorithm is EEkt, power `p` must be specified.
If the algorithm is Valencia, you can provide `p` (β) and `γ`, or pass `β` explicitly to override `p`.
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

    # Compute invR2 once and thread it through
    invR2 = inv(R * R)
    # Now call the unified implementation with conditional logic.
    return _ee_genkt_algorithm(particles = recombination_particles, p = p, R = R,
                               invR2 = invR2, algorithm = algorithm, recombine = recombine,
                               γ = γ)
end

"""
    _ee_genkt_algorithm(particles::AbstractVector{EEJet};
                        algorithm::JetAlgorithm.Algorithm, p::Real, R::Real,
                        invR2::Union{Real, Nothing} = nothing, recombine = addjets, γ::Real = 1.0,
                        beta::Union{Real, Nothing} = nothing)

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
 - `γ::Real = 1.0`: Angular exponent for the Valencia beam metric (ignored for other algorithms).
 - `beta::Union{Real, Nothing} = nothing`: Optional alias for the Valencia energy exponent (β).
     When provided with `algorithm == JetAlgorithm.Valencia`, it overrides `p`.

# Returns
- `clusterseq`: The resulting `ClusterSequence` object representing the
  reconstructed jets.
"""
function _ee_genkt_algorithm(; particles::AbstractVector{EEJet},
                             algorithm::JetAlgorithm.Algorithm, p::Real, R::Real,
                             invR2::Union{Real, Nothing} = nothing, recombine = addjets,
                             γ::Real = 1.0,
                             beta::Union{Real, Nothing} = nothing)
    # Bounds
    N::Int = length(particles)

    # invR2 provided by caller when available; otherwise compute from R once here
    if invR2 === nothing
        invR2 = inv(R * R)
    end
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
    elseif algorithm == JetAlgorithm.Valencia
        dij_factor = 1.0
    else
        throw(ArgumentError("Algorithm $algorithm not supported for e+e-"))
    end

    # For optimised reconstruction generate an SoA containing the necessary
    # jet information and populate it accordingly
    # We need N slots for this array
    eereco = StructArray{EERecoJet}(undef, N)

    fill_reco_array!(eereco, particles, invR2, p)

    # Setup the initial history and get the total energy
    history, Qtot = initial_history(particles)

    clusterseq = ClusterSequence(algorithm, p, R, RecoStrategy.N2Plain, particles, history,
                                 Qtot)

    # Run over initial pairs of jets to find nearest neighbours
    get_angular_nearest_neighbours!(eereco, algorithm, dij_factor, p, invR2, γ)

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
            insert_new_jet!(eereco, ijetA, newjet_k, invR2, merged_jet, p)
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
                    update_nn_no_cross!(eereco, i, N, algorithm, dij_factor, invR2, p, γ)
                end
            end
        end

        # Finally, we need to update the nearest neighbours for the new jet, checking both ways
        # (But only if there was a new jet!)
        if ijetA != ijetB
            update_nn_cross!(eereco, ijetA, N, algorithm, dij_factor, invR2, p, γ)
        end

        # Only for debugging purposes...
        # ee_check_consistency(clusterseq, eereco, N)
    end

    # Return the final cluster sequence structure
    clusterseq
end
