using LoopVectorization

"""
    dist(i, j, rapidity_array, phi_array)

Compute the distance between points in a 2D space defined by rapidity and phi coordinates.

# Arguments
- `i::Int`: Index of the first point to consider (indexes into `rapidity_array` and `phi_array`).
- `j::Int`: Index of the second point to consider (indexes into `rapidity_array` and `phi_array`).
- `rapidity_array::Vector{Float64}`: Array of rapidity coordinates.
- `phi_array::Vector{Float64}`: Array of phi coordinates.

# Returns
- `distance::Float64`: The distance between the two points.
"""
Base.@propagate_inbounds function dist(i, j, rapidity_array, phi_array)
    drapidity = rapidity_array[i] - rapidity_array[j]
    dphi = abs(phi_array[i] - phi_array[j])
    dphi = ifelse(dphi > pi, 2pi - dphi, dphi)
    @muladd drapidity * drapidity + dphi * dphi
end

"""
    dij(i, kt2_array, nn, nndist)

Compute the dij value for a given index `i` to its nearest neighbor. The nearest
neighbor is determined from `nn[i]`, and the metric distance to the nearest
neighbor is given by the distance `nndist[i]` applying the lower of the
`kt2_array` values for the two particles.ßß

# Arguments
- `i`: The index of the element.
- `kt2_array`: An array of kt2 values.
- `nn`: An array of nearest neighbors.
- `nndist`: An array of nearest neighbor distances.

# Returns
- The computed dij value.
"""
Base.@propagate_inbounds function dij(i, kt2_array, nn, nndist)
    j = nn[i]
    d = nndist[i]
    d * min(kt2_array[i], kt2_array[j])
end

"""
    upd_nn_crosscheck!(i, from, to, rapidity_array, phi_array, R2, nndist, nn)

Update the nearest neighbor information for a given particle index `i` against
all particles in the range indexes `from` to `to`. The function updates the
`nndist` and `nn` arrays with the nearest neighbor distance and index
respectively, both for particle `i` and the checked particles `[from:to]` (hence
*crosscheck*).

# Arguments
- `i::Int`: The index of the particle to update and check against.
- `from::Int`: The starting index of the range of particles to check against.
- `to::Int`: The ending index of the range of particles to check against.
- `rapidity_array`: An array containing the rapidity values of all particles.
- `phi_array`: An array containing the phi values of the all particles.
- `R2`: The squared jet distance threshold for considering a particle as a
  neighbour.
- `nndist`: The array that stores the nearest neighbor distances.
- `nn`: The array that stores the nearest neighbor indices.
"""
Base.@propagate_inbounds function upd_nn_crosscheck!(i::Int, from::Int, to::Int,
                                                     rapidity_array, phi_array, R2, nndist,
                                                     nn)
    nndist_min = R2
    nn_min = i
    @inbounds @simd for j in from:to
        Δ2 = dist(i, j, rapidity_array, phi_array)
        if Δ2 < nndist_min
            nn_min = j
            nndist_min = Δ2
        end
        if Δ2 < nndist[j]
            nndist[j] = Δ2
            nn[j] = i
        end
    end
    nndist[i] = nndist_min
    nn[i] = nn_min
end

# finds new nn for i

"""
    upd_nn_nocross!(i, from, to, rapidity_array, phi_array, R2, nndist, nn)

Update the nearest neighbor information for a given particle index `i` against
all particles in the range indexes `from` to `to`. The function updates the
`nndist` and `nn` arrays with the nearest neighbor distance and index
respectively, only for particle `i` (hence *nocross*).

# Arguments
- `i::Int`: The index of the particle to update and check against.
- `from::Int`: The starting index of the range of particles to check against.
- `to::Int`: The ending index of the range of particles to check against.
- `rapidity_array`: An array containing the rapidity values of all particles.
- `phi_array`: An array containing the phi values of the all particles.
- `R2`: The squared jet distance threshold for considering a particle as a
  neighbour.
- `nndist`: The array that stores the nearest neighbor distances.
- `nn`: The array that stores the nearest neighbor indices.
"""
Base.@propagate_inbounds function upd_nn_nocross!(i::Int, from::Int, to::Int,
                                                  rapidity_array, phi_array, R2, nndist, nn)
    nndist_min = R2
    nn_min = i
    @inbounds @simd for j in from:(i - 1)
        Δ2 = dist(i, j, rapidity_array, phi_array)
        if Δ2 <= nndist_min
            nn_min = j
            nndist_min = Δ2
        end
    end
    @inbounds @simd for j in (i + 1):to
        Δ2 = dist(i, j, rapidity_array, phi_array)
        f = Δ2 <= nndist_min
        nn_min = ifelse(f, j, nn_min)
        nndist_min = ifelse(f, Δ2, nndist_min)
    end
    nndist[i] = nndist_min
    nn[i] = nn_min
end

"""
    upd_nn_step!(i, j, k, N, Nn, kt2_array, rapidity_array, phi_array, R2, nndist, nn, nndij)

Update the nearest neighbor information after a jet merge step.

Arguments:
- `i`: Index of the first particle in the last merge step.
- `j`: Index of the second particle in the last merge step.
- `k`: Index of the current particle for which the nearest neighbour will be
  updated.
- `N`: Total number of particles (currently vaild array indexes are `[1:N]`).
- `Nn`: Number of nearest neighbors to consider.
- `kt2_array`: Array of transverse momentum squared values.
- `rapidity_array`: Array of rapidity values.
- `phi_array`: Array of azimuthal angle values.
- `R2`: Distance threshold squared for nearest neighbors.
- `nndist`: Array of nearest neighbor geometric distances.
- `nn`: Array of nearest neighbor indices.
- `nndij`: Array of metric distances between particles.

This function updates the nearest neighbor information for the current particle
`k` by considering the distances to particles `i` and `j`. It checks if the
distance between `k` and `i` is smaller than the current nearest neighbor
distance for `k`, and updates the nearest neighbor information accordingly. It
also updates the nearest neighbor information for `i` if the distance between
`k` and `i` is smaller than the current nearest neighbor distance for `i`.
Finally, it checks if the nearest neighbor of `k` is the total number of
particles `Nn` and updates it to `j` if necessary.

"""
Base.@propagate_inbounds function upd_nn_step!(i, j, k, N, Nn, kt2_array, rapidity_array,
                                               phi_array, R2, nndist, nn, nndij)
    nnk = nn[k] # Nearest neighbour of k
    if nnk == i || nnk == j
        # Our old nearest neighbour is one of the merged particles
        upd_nn_nocross!(k, 1, N, rapidity_array, phi_array, R2, nndist, nn) # Update dist and nn
        nndij[k] = dij(k, kt2_array, nn, nndist)
        nnk = nn[k]
    end

    if j != i && k != i
        # Update particle's nearest neighbour if it's not i and the merge step was not a beam merge
        Δ2 = dist(i, k, rapidity_array, phi_array)
        if Δ2 < nndist[k]
            nndist[k] = Δ2
            nnk = nn[k] = i
            nndij[k] = dij(k, kt2_array, nn, nndist)
        end

        cond = Δ2 < nndist[i]
        nndist[i], nn[i] = ifelse(cond, (Δ2, k), (nndist[i], nn[i]))
    end

    # If the previous nearest neighbour was the final jet in the array before
    # the merge that was just done, this jet has now been moved in the array to
    # position k (to compactify the array), so we need to update the nearest
    # neighbour
    nnk == Nn && (nn[k] = j)
end

"""
    plain_jet_reconstruct(particles::Vector{T}; p = -1, R = 1.0, recombine = +) where T

Perform pp jet reconstruction using the plain algorithm.

# Arguments
- `particles::Vector{T}`: A vector of particles used for jet reconstruction, any
   array of particles, which supports suitable 4-vector methods, viz. pt2(),
   phi(), rapidity(), px(), py(), pz(), energy(), can be used. for each element.
- `algorithm::Union{JetAlgorithm, Nothing} = nothing`: The explicit jet
  algorithm to use.
- `p::Int=-1`: The integer value used for jet reconstruction.
- `R::Float64=1.0`: The radius parameter used for jet reconstruction.
- `recombine::Function=+`: The recombination function used for jet
  reconstruction.

**Note** for the `particles` argument, the 4-vector methods need to exist in the
JetReconstruction package namespace.

This code will use the `k_t` algorithm types, operating in `(rapidity, φ)` space.

It is not necessary to specify both the `algorithm` and the `p` (power) value.
If both are given they must be consistent or an exception is thrown.

# Returns
- `Vector{PseudoJet}`: A vector of reconstructed jets.

# Example
```julia
jets = plain_jet_reconstruct(particles; p = -1, R = 0.4)
jets = plain_jet_reconstruct(particles; algorithm = JetAlgorithm.Kt, R = 1.0)
```
"""
function plain_jet_reconstruct(particles::Vector{T}; p::Union{Real, Nothing} = -1, R = 1.0,
                               algorithm::Union{JetAlgorithm.Algorithm, Nothing} = nothing,
                               recombine = +) where {T}

    # Check for consistency between algorithm and power
    (p, algorithm) = get_algorithm_power_consistency(p = p, algorithm = algorithm)

    # Integer p if possible
    p = (round(p) == p) ? Int(p) : p

    if T == PseudoJet
        # recombination_particles will become part of the cluster sequence, so size it for
        # the starting particles and all N recombinations
        recombination_particles = copy(particles)
        sizehint!(recombination_particles, length(particles) * 2)
    else
        recombination_particles = PseudoJet[]
        sizehint!(recombination_particles, length(particles) * 2)
        for i in eachindex(particles)
            push!(recombination_particles,
                  PseudoJet(px(particles[i]), py(particles[i]), pz(particles[i]),
                            energy(particles[i])))
        end
    end

    # Now call the actual reconstruction method, tuned for our internal EDM
    _plain_jet_reconstruct(particles = recombination_particles, p = p, R = R,
                           algorithm = algorithm,
                           recombine = recombine)
end

"""
    _plain_jet_reconstruct(; particles::Vector{PseudoJet}, p = -1, R = 1.0, recombine = +)

This is the internal implementation of jet reconstruction using the plain
algorithm. It takes a vector of `particles` representing the input particles and
reconstructs jets based on the specified parameters. Here the particles must be
of type `PseudoJet`.

Users of the package should use the `plain_jet_reconstruct` function as their
entry point to this jet reconstruction.

The power value maps to specific pp jet reconstruction algorithms: -1 = AntiKt,
0 = Cambridge/Aachen, 1 = Inclusive Kt. Floating point values are allowed for
generalised k_t algorithm.

# Arguments
- `particles`: A vector of `PseudoJet` objects representing the input particles.
- `p=-1`: The power to which the transverse momentum (`pt`) of each particle is
  raised.
- `R=1.0`: The jet radius parameter.
- `recombine`: The recombination function used to merge two jets. Default is `+`
  (additive recombination).

# Returns
- `clusterseq`: The resulting `ClusterSequence` object representing the
  reconstructed jets.
"""
function _plain_jet_reconstruct(; particles::Vector{PseudoJet}, p = -1, R = 1.0,
                                algorithm::JetAlgorithm.Algorithm = JetAlgorithm.AntiKt,
                                recombine = +)
    # Bounds
    N::Int = length(particles)
    # Parameters
    R2 = R^2

    # Optimised compact arrays for determining the next merge step
    # We make sure these arrays are type stable - have seen issues where, depending on the values
    # returned by the methods, they can become unstable and performance degrades
    kt2_array::Vector{Float64} = pt2.(particles) .^ p
    phi_array::Vector{Float64} = phi.(particles)
    rapidity_array::Vector{Float64} = rapidity.(particles)
    nn::Vector{Int} = Vector(1:N) # nearest neighbours
    nndist::Vector{Float64} = fill(float(R2), N) # geometric distances to the nearest neighbour
    nndij::Vector{Float64} = zeros(N) # dij metric distance

    # Maps index from the compact array to the clusterseq jet vector
    clusterseq_index::Vector{Int} = collect(1:N)

    # Setup the initial history and get the total energy
    history, Qtot = initial_history(particles)
    # Current implementation mutates the particles vector, so need to copy it
    # for the cluster sequence (there is too much copying happening, so this
    # needs to be rethought and reoptimised)
    clusterseq = ClusterSequence(algorithm, p, RecoStrategy.N2Plain, particles, history,
                                 Qtot)

    # Initialize nearest neighbours
    @simd for i in 1:N
        upd_nn_crosscheck!(i, 1, i - 1, rapidity_array, phi_array, R2, nndist, nn)
    end

    # diJ table * R2
    @inbounds @simd for i in 1:N
        nndij[i] = dij(i, kt2_array, nn, nndist)
    end

    iteration::Int = 1
    while N != 0
        # Extremely odd - having these @debug statements present causes a performance
        # degradation of ~140μs per event on my M2 mac (20%!), even when no debugging is used
        # so they need to be completely commented out...
        #@debug "Beginning iteration $iteration"

        # Findmin and add back renormalisation to distance
        dij_min, i = fast_findmin(nndij, N)
        @fastmath dij_min /= R2
        j::Int = nn[i]

        #@debug "Closest compact jets are $i ($(clusterseq_index[i])) and $j ($(clusterseq_index[j]))"

        if i != j # Merge jets i and j
            # swap if needed
            if j < i
                i, j = j, i
            end

            # Source "history" for merge
            hist_i = clusterseq.jets[clusterseq_index[i]]._cluster_hist_index
            hist_j = clusterseq.jets[clusterseq_index[j]]._cluster_hist_index

            # Recombine i and j into the next jet
            push!(clusterseq.jets,
                  recombine(clusterseq.jets[clusterseq_index[i]],
                            clusterseq.jets[clusterseq_index[j]]))
            # Get its index and the history index
            newjet_k = length(clusterseq.jets)
            newstep_k = length(clusterseq.history) + 1
            clusterseq.jets[newjet_k]._cluster_hist_index = newstep_k
            # Update history
            add_step_to_history!(clusterseq, minmax(hist_i, hist_j)..., newjet_k, dij_min)

            # Update the compact arrays, reusing the i-th slot
            kt2_array[i] = pt2(clusterseq.jets[newjet_k])^p
            rapidity_array[i] = rapidity(clusterseq.jets[newjet_k])
            phi_array[i] = phi(clusterseq.jets[newjet_k])
            clusterseq_index[i] = newjet_k
            nndist[i] = R2
            nn[i] = i
        else # i == j, this is a final jet ("merged with beam")
            add_step_to_history!(clusterseq,
                                 clusterseq.jets[clusterseq_index[i]]._cluster_hist_index,
                                 BeamJet, Invalid, dij_min)
        end

        # Squash step - copy the final jet's compact data into the j-th slot
        if j != N
            phi_array[j] = phi_array[N]
            rapidity_array[j] = rapidity_array[N]
            kt2_array[j] = kt2_array[N]
            nndist[j] = nndist[N]
            nn[j] = nn[N]
            nndij[j] = nndij[N]
            clusterseq_index[j] = clusterseq_index[N]
        end

        Nn::Int = N
        N -= 1
        iteration += 1

        # Update nearest neighbours step
        @inbounds @simd for k in 1:N
            upd_nn_step!(i, j, k, N, Nn, kt2_array, rapidity_array, phi_array, R2, nndist,
                         nn, nndij)
        end

        nndij[i] = dij(i, kt2_array, nn, nndist)
    end

    # Return the final cluster sequence structure
    clusterseq
end
