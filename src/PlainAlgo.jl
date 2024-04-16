using LoopVectorization

Base.@propagate_inbounds function dist(i, j, rapidity_array, phi_array)
    drapidity = rapidity_array[i] - rapidity_array[j]
    dphi = abs(phi_array[i] - phi_array[j])
    dphi = ifelse(dphi > pi, 2pi - dphi, dphi)
    muladd(drapidity, drapidity, dphi * dphi)
end

# d_{ij} distance with i's NN (times R^2)
Base.@propagate_inbounds function dij(i, kt2_array, nn, nndist)
    j = nn[i]
    d = nndist[i]
    d * min(kt2_array[i], kt2_array[j])
end

# finds new nn for i and checks everyone additionally
Base.@propagate_inbounds function upd_nn_crosscheck!(i::Int, from::Int, to::Int, rapidity_array, phi_array, R2, nndist, nn)
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
Base.@propagate_inbounds function upd_nn_nocross!(i::Int, from::Int, to::Int, rapidity_array, phi_array, R2, nndist, nn)
    nndist_min = R2
    nn_min = i
    @inbounds @simd for j in from:(i-1)
        Δ2 = dist(i, j, rapidity_array, phi_array)
        if Δ2 <= nndist_min
            nn_min = j
            nndist_min = Δ2
        end
    end
    @inbounds @simd for j in (i+1):to
        Δ2 = dist(i, j, rapidity_array, phi_array)
        f = Δ2 <= nndist_min
        nn_min = ifelse(f, j, nn_min)
        nndist_min = ifelse(f, Δ2, nndist_min)
    end
    nndist[i] = nndist_min
    nn[i] = nn_min
end

# entire NN update step
Base.@propagate_inbounds function upd_nn_step!(i, j, k, N, Nn, kt2_array, rapidity_array, phi_array, R2, nndist, nn, nndij)
    nnk = nn[k]
    if nnk == i || nnk == j
        upd_nn_nocross!(k, 1, N, rapidity_array, phi_array, R2, nndist, nn) # update dist and nn
        nndij[k] = dij(k, kt2_array, nn, nndist)
        nnk = nn[k]
    end

    if j != i && k != i
        Δ2 = dist(i, k, rapidity_array, phi_array)
        if Δ2 < nndist[k]
            nndist[k] = Δ2
            nnk = nn[k] = i
            nndij[k] = dij(k, kt2_array, nn, nndist)
        end

        cond = Δ2 < nndist[i]
        nndist[i], nn[i] = ifelse(cond, (Δ2, k), (nndist[i], nn[i]))
    end

    nnk == Nn && (nn[k] = j)
end

"""
This is the N2Plain jet reconstruction algorithm interface, called with an arbitrary array
of particles, which supports suitable 4-vector methods, viz.
 - pt2(), phi(), rapidity(), px(), py(), pz(), energy()
for each element.

Examples of suitable types are JetReconstruction.PseudoJet and LorentzVectorHEP.
N.B. these methods need to exist in the namespace of this package, i.e. JetReconstruction.pt2(::T),
which is already done for the two types above.
"""
function plain_jet_reconstruct(particles::Vector{T}; p = -1, R = 1.0, recombine = +) where T
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
            push!(recombination_particles, PseudoJet(px(particles[i]), py(particles[i]), pz(particles[i]), energy(particles[i])))
        end
    end

    # Now call the actual reconstruction method, tuned for our internal EDM
    _plain_jet_reconstruct(particles=recombination_particles, p=p, R=R, recombine=recombine)
end


function _plain_jet_reconstruct(;particles::Vector{PseudoJet}, p = -1, R = 1.0, recombine = +)
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
    nn = Vector(1:N) # nearest neighbours
    nndist = fill(float(R2), N) # geometric distances to the nearest neighbour
    nndij::Vector{Float64} = zeros(N) # dij metric distance

    # Maps index from the compact array to the clusterseq jet vector
    clusterseq_index::Vector{Int} = collect(1:N)

    # Setup the initial history and get the total energy
    history, Qtot = initial_history(particles)
    # Current implementation mutates the particles vector, so need to copy it
    # for the cluster sequence (there is too much copying happening, so this
    # needs to be rethought and reoptimised)
    clusterseq = ClusterSequence(particles, history, Qtot)

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
        @debug "Beginning iteration $iteration"
        # Findmin and add back renormalisation to distance
        dij_min, i = fast_findmin(nndij, N)
        dij_min *= R2
        j::Int = nn[i]
        @debug "Closest compact jets are $i ($(clusterseq_index[i])) and $j ($(clusterseq_index[j]))"

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
                recombine(clusterseq.jets[clusterseq_index[i]], clusterseq.jets[clusterseq_index[j]]))
            # Get its index and the history index
            newjet_k = length(clusterseq.jets)
            newstep_k = length(clusterseq.history) + 1
            clusterseq.jets[newjet_k]._cluster_hist_index = newstep_k
            # Update history
            add_step_to_history!(clusterseq, minmax(hist_i, hist_j)..., newjet_k, dij_min)

            # Update the compact arrays, reusing the i-th slot
            kt2_array[i] = pt2(clusterseq.jets[newjet_k]) ^ p
            rapidity_array[i] = rapidity(clusterseq.jets[newjet_k])
            phi_array[i] = phi(clusterseq.jets[newjet_k])
            clusterseq_index[i] = newjet_k
            nndist[i] = R2
            nn[i] = i
        else # i == j, this is a final jet ("merged with beam")
            add_step_to_history!(clusterseq, clusterseq.jets[clusterseq_index[i]]._cluster_hist_index, BeamJet, Invalid, dij_min)
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
            upd_nn_step!(i, j, k, N, Nn, kt2_array, rapidity_array, phi_array, R2, nndist, nn, nndij)
        end

        nndij[i] = dij(i, kt2_array, nn, nndist)
    end

    # Return the final cluster sequence structure
    clusterseq
end
