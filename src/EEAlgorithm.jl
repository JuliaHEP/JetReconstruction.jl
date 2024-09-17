# Use these constants whenever we need to set a large value for the distance or
# metric distance
const large_distance = 16.0 # = 4^2
const large_dij = 1.0e6

"""
    angular_distance(jet1::EEjet, jet2::EEjet) -> Float64

Calculate the angular distance between two jets using the formula (1 - cos(θ)).

# Arguments
- `jet1::EEjet`: The first jet.
- `jet2::EEjet`: The second jet.

# Returns
- `Float64`: The angular distance between `jet1` and `jet2`, which is ``1 - cos
  \theta``.
"""
@inline function angular_distance(jet1::EEjet, jet2::EEjet)
    # Calculate the angular distance between two jets (1 - cos(θ))
    @muladd 1.0 - nx(jet1) * nx(jet2) - ny(jet1) * ny(jet2) - nz(jet1) * nz(jet2)
end

"""
    dij_dist(nndist, jet1::EEjet, jet2::EEjet, p = 1)

Calculate the dij distance between two ``e^+e^-``jets.

# Arguments
- `nndist`: The angular nearest neighbor distance.
- `jet1::EEjet`: The first jet object of type `EEjet`.
- `jet2::EEjet`: The second jet object of type `EEjet`.
- `p=1`: The power weighting for the energy of the jets.

# Returns
- The dij distance between `jet1` and `jet2`.

# Notes
- The function uses the energy of the jets raised to the power of `2p` and takes the minimum of these values multiplied by `nndist`.
"""
@inline function dij_dist(nndist, jet1::EEjet, jet2::EEjet, p = 1)
    # Calculate the dij distance between two jets
    nndist * min(energy(jet1)^2p, energy(jet2)^2p)
end

function get_angular_nearest_neighbours!(cs::ClusterSequence,
                                         clusterseq_index::Vector{Int},
                                         nndist::Vector{Float64},
                                         nndij::Vector{Float64}, nni::Vector{Int})
    jets = cs.jets
    # Get the nearest neighbour for each jet
    @inbounds for i in eachindex(jets)
        @inbounds for j in (i + 1):length(jets)
            this_nndist = angular_distance(jets[clusterseq_index[i]],
                                           jets[clusterseq_index[j]])
            if this_nndist < nndist[i]
                nndist[i] = this_nndist
                nni[i] = j
            end
            if this_nndist < nndist[j]
                nndist[j] = this_nndist
                nni[j] = i
            end
        end
    end
    @inbounds for i in eachindex(jets)
        nndij[i] = dij_dist(nndist[i], jets[clusterseq_index[i]],
                            jets[clusterseq_index[nni[i]]], cs.power)
    end
end

function get_angular_nearest_neighbours!(eereco, algorithm, dij_factor)
    # Get the initial nearest neighbours for each jet
    N = length(eereco)
    this_dist_vector = Vector{Float64}(undef, N)
    # Nearest neighbour geometric distance
    @inbounds for i in 1:N
        # TODO: Replace the 'j' loop with a vectorised operation over the appropriate array elements
        # this_dist_vector .= 1.0 .- eereco.nx[i:N] .* eereco[i + 1:end].nx .-
        #     eereco[i].ny .* eereco[i + 1:end].ny .- eereco[i].nz .* eereco[i + 1:end].nz
        @inbounds for j in (i + 1):N
            @muladd this_nndist = 1.0 - eereco[i].nx * eereco[j].nx - eereco[i].ny * eereco[j].ny - eereco[i].nz * eereco[j].nz
            if this_nndist < eereco[i].nndist
                eereco[i].nndist = this_nndist
                eereco[i].nni = j
            end
            if this_nndist < eereco[j].nndist
                eereco[j].nndist = this_nndist
                eereco[j].nni = i
            end
        end
    end
    # Nearest neighbour dij distance
    @inbounds for i in 1:N
        eereco[i].dijdist = min(eereco[i].E2p, eereco[eereco[i].nni].E2p) * dij_factor * eereco[i].nndist
    end
    # For the EEKt algorithm, we need to check the beam distance as well
    # This is structured to only check for EEKt once!
    if algorithm == JetAlgorithm.EEKt
        @inbounds for i in 1:N
            if eereco[i].E2p < eereco[i].dijdist
                eereco[i].dijdist = eereco[i].E2p
                eereco[i].nni = 0
            end
        end
    end
end

function update_nn_no_cross!(i, N, cs::ClusterSequence, clusterseq_index, nndist, nndij,
                             nni)
    # Update the nearest neighbour for jet i, w.r.t. all other active jets
    nndist[i] = large_distance
    nni[i] = i
    jets = cs.jets
    @inbounds for j in 1:N
        if j != i
            this_nndist = angular_distance(jets[clusterseq_index[i]],
                                           jets[clusterseq_index[j]])
            if this_nndist < nndist[i]
                nndist[i] = this_nndist
                nni[i] = j
            end
        end
    end
    nndij[i] = dij_dist(nndist[i], jets[clusterseq_index[i]],
                        jets[clusterseq_index[nni[i]]], cs.power)
end

# Update the nearest neighbour for jet i, w.r.t. all other active jets
function update_nn_no_cross!(eereco, i, N, algorithm, dij_factor)
    eereco[i].nndist = large_distance
    eereco[i].nni = 0
    @inbounds for j in 1:N
        if j != i
            @muladd this_nndist = 1.0 - eereco[i].nx * eereco[j].nx - eereco[i].ny * eereco[j].ny - eereco[i].nz * eereco[j].nz
            if this_nndist < eereco[i].nndist
                eereco[i].nndist = this_nndist
                eereco[i].nni = j
            end
        end
    end
    eereco[i].dijdist = min(eereco[i].E2p, eereco[eereco[i].nni].E2p) * dij_factor * eereco[i].nndist
    if algorithm == JetAlgorithm.EEKt
        if eereco[i].E2p < eereco[i].dijdist
            eereco[i].dijdist = eereco[i].E2p
            eereco[i].nni = 0
        end
    end
end


function update_nn_cross!(i, N, cs::ClusterSequence, clusterseq_index, nndist, nndij, nni)
    # Update the nearest neighbour for jet i, w.r.t. all other active jets
    # also doing the cross check for the other jet
    nndist[i] = large_distance
    nni[i] = i
    jets = cs.jets
    @inbounds for j in 1:N
        if j != i
            this_nndist = angular_distance(jets[clusterseq_index[i]],
                                           jets[clusterseq_index[j]])
            if this_nndist < nndist[i]
                nndist[i] = this_nndist
                nni[i] = j
            end
            if this_nndist < nndist[j]
                nndist[j] = this_nndist
                nni[j] = i
                # j will not be revisited, so update metric distance here
                nndij[j] = dij_dist(nndist[j], jets[clusterseq_index[j]],
                                    jets[clusterseq_index[i]], cs.power)
            end
        end
    end
    nndij[i] = dij_dist(nndist[i], jets[clusterseq_index[i]],
                        jets[clusterseq_index[nni[i]]], cs.power)
end

function update_nn_cross!(eereco, i, N, algorithm, dij_factor)
    # Update the nearest neighbour for jet i, w.r.t. all other active jets
    # also doing the cross check for the other jet
    eereco[i].nndist = large_distance
    eereco[i].nni = 0
    @inbounds for j in 1:N
        if j != i
            @muladd this_nndist = 1.0 - eereco[i].nx * eereco[j].nx - eereco[i].ny * eereco[j].ny - eereco[i].nz * eereco[j].nz
            if this_nndist < eereco[i].nndist
                eereco[i].nndist = this_nndist
                eereco[i].nni = j
            end
            if this_nndist < eereco[j].nndist
                eereco[j].nndist = this_nndist
                eereco[j].nni = i
                # j will not be revisited, so update metric distance here
                nndij[j] = min(eereco[i].E2p, eereco[j].E2p) * dij_factor * eereco[j].nndist
                if algorithm == JetAlgorithm.EEKt
                    if eereco[j].E2p < eereco[j].dijdist
                        eereco[j].dijdist = eereco[j].E2p
                        eereco[j].nni = 0
                    end
                end
            end
        end
    end
    eereco[i].dijdist = min(eereco[i].E2p, eereco[eereco[i].nni].E2p) * dij_factor * eereco[i].nndist
    if algorithm == JetAlgorithm.EEKt
        if eereco[i].E2p < eereco[i].dijdist
            eereco[i].dijdist = eereco[i].E2p
            eereco[i].nni = 0
        end
    end
end

function ee_check_consistency(clusterseq, clusterseq_index, N, nndist, nndij, nni, msg)
    # Check the consistency of the reconstruction state
    for i in 1:N
        if nni[i] > N
            @error "Jet $i has invalid nearest neighbour $(nni[i])"
        end
        for i in 1:N
            jet = clusterseq.jets[clusterseq_index[i]]
            jet_hist = clusterseq.history[jet._cluster_hist_index]
            if jet_hist.child != Invalid
                @error "Jet $i has invalid child $(jet_hist.child)"
            end
        end
    end
    @debug "Consistency check passed at $msg"
end

function dij_correct_for_beam!(cs::ClusterSequence, clusterseq_index, nndist, nndij, nni, N,
                               dij_factor)
    # Correct the dij metric for the beam merges
    for i in 1:N
        diB = energy(cs.jets[clusterseq_index[i]])^(2 * cs.power)
        if diB < nndij[i] * dij_factor
            nndij[i] = diB / dij_factor
            nni[i] = 0
        end
    end
end

function ee_genkt_algorithm(particles::Vector{T}; p::Union{Real, Nothing} = -1, R = 4.0,
                            algorithm::Union{JetAlgorithm.Algorithm, Nothing} = nothing,
                            recombine = +) where {T}

    # Check for consistency between algorithm and power
    (p, algorithm) = get_algorithm_power_consistency(p = p, algorithm = algorithm)

    # Integer p if possible
    p = (round(p) == p) ? Int(p) : p

    # For the Durham algorithm, p=1 and R is not used, but nominally set to 4
    if algorithm == JetAlgorithm.Durham
        R = 4.0
    end

    if T == EEjet
        # recombination_particles will become part of the cluster sequence, so size it for
        # the starting particles and all N recombinations
        recombination_particles = copy(particles)
        sizehint!(recombination_particles, length(particles) * 2)
    else
        recombination_particles = EEjet[]
        sizehint!(recombination_particles, length(particles) * 2)
        for i in eachindex(particles)
            push!(recombination_particles,
                  EEjet(px(particles[i]), py(particles[i]), pz(particles[i]),
                        energy(particles[i])))
        end
    end

    # Now call the actual reconstruction method, tuned for our internal EDM
    _ee_genkt_algorithm(particles = recombination_particles, p = p, R = R,
                        algorithm = algorithm,
                        recombine = recombine)
end

function _ee_genkt_algorithm(; particles::Vector{EEjet}, p = 1, R = 4.0,
                             algorithm::JetAlgorithm.Algorithm = JetAlgorithm.Durham,
                             recombine = +)
    # Bounds
    N::Int = length(particles)

    # R squared
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
    for i in eachindex(particles)
        eereco.index[i] = i
        eereco.nni[i] = 0
        eereco.nndist[i] = R2
        # eereco[i].dijdist = UNDEF # Not needed
        eereco.nx[i] = nx(particles[i])
        eereco.ny[i] = ny(particles[i])
        eereco.nz[i] = nz(particles[i])
        eereco.E2p[i] = energy(particles[i])^(2p)
    end
    # println("Jet 23 compact init: ", eereco[23])
    # println("Jet 23 particles: ", particles[23])

    # Optimised compact arrays for determining the next merge step We make sure
    # these arrays are type stable - have seen issues where, depending on the
    # values returned by the methods, they can become unstable and performance
    # degrades
    # nndist::Vector{Float64} = Vector{Float64}(undef, N) # distances 
    # fill!(nndist, R2)
    # nndij::Vector{Float64} = Vector{Float64}(undef, N)  # dij metric distance
    # nni::Vector{Int} = collect(1:N) # Nearest neighbour index (in the compact arrays!)

    # # Maps index from the compact array to the clusterseq jet vector
    # clusterseq_index::Vector{Int} = collect(1:N)

    # Setup the initial history and get the total energy
    history, Qtot = initial_history(particles)

    clusterseq = ClusterSequence(algorithm, p, R, RecoStrategy.N2Plain, particles, history,
                                 Qtot)
    # println("Jet 23: ", clusterseq.jets[23])

    # Run over initial pairs of jets to find nearest neighbours
    get_angular_nearest_neighbours!(eereco, algorithm, dij_factor)

    # println("Jet 23 compact: ", eereco[23])

    return clusterseq

    # ee_check_consistency(clusterseq, clusterseq_index, N, nndist, nndij, nni, "Start")

    # Now we can start the main loop
    iter = 0
    while N != 0
        iter += 1

        # if algorithm == JetAlgorithm.EEKt
        #     dij_correct_for_beam!(clusterseq, clusterseq_index, nndist, nndij, nni, N,
        #                           dij_factor)
        # end

        dij_min, ijetA = fast_findmin(eereco.dijdist, N)
        ijetB = eereco[ijetA].nni

        println("N $N, Iteration $iter: dij_min = $dij_min, ijetA = $ijetA, ijetB = $ijetB; njets = $(length(clusterseq.jets))")

        # Normalise the dij_min
        # dij_min *= dij_factor

        # Now we check if there is a "beam" merge possibility
        if ijetB == 0
            # We have an EEKt beam merge
            ijetB == ijetA
            println("Beam merge: dij_min = $dij_min, ijetA = $ijetA, ijetB = $ijetB")
            add_step_to_history!(clusterseq,
                                 clusterseq.jets[eereco[ijetA].index]._cluster_hist_index,
                                 BeamJet, Invalid, dij_min)
        elseif N == 1
            # We have a single jet left
            add_step_to_history!(clusterseq,
                                 clusterseq.jets[eereco[ijetA].index]._cluster_hist_index,
                                 BeamJet, Invalid, dij_min)
        else
            # Jet-jet merge
            if ijetB < ijetA
                ijetA, ijetB = ijetB, ijetA
            end

            # Source "history" for merge
            hist_jetA = clusterseq.jets[eereco[ijetA].index]._cluster_hist_index
            hist_jetB = clusterseq.jets[eereco[ijetB].index]._cluster_hist_index

            # Recombine jetA and jetB into the next jet
            merged_jet = recombine(clusterseq.jets[eereco[ijetA].index],
                                   clusterseq.jets[eereco[ijetB].index])
            merged_jet._cluster_hist_index = length(clusterseq.history) + 1

            # Now add the jet to the sequence, and update the history
            push!(clusterseq.jets, merged_jet)
            newjet_k = length(clusterseq.jets)
            add_step_to_history!(clusterseq, minmax(hist_jetA, hist_jetB)...,
                                 newjet_k, dij_min)

            # Update the compact arrays, reusing the JetA slot
            eereco[ijetA].index = newjet_k
            eereco[ijetA].nni = 0
            eereco[ijetA].nndist = R2
            eereco[ijetA].nx = nx(merged_jet)
            eereco[ijetA].ny = ny(merged_jet)
            eereco[ijetA].nz = nz(merged_jet)
            eereco[ijetA].E2p = energy(merged_jet)^(2p)
        end

        # Squash step - copy the final jet's compact data into the jetB slot
        if ijetB != N
            eereco[ijetB] = eereco[N]
        end

        # Now number of active jets is decreased by one
        N -= 1

        # Update nearest neighbours step
        for i in 1:N
            # First, if any jet's NN was the old index N, it's now ijetB
            if (ijetB != N + 1) && (eereco[i].nni == N + 1)
                eereco[i].nni = ijetB
            else
                # Otherwise, if the jet had ijetA or ijetB as their NN, we need to update them
                # plus "belt and braces" check for an invalid NN (>N)
                if (eereco[i].nni == ijetA) || (eereco[i].nni == ijetB) || (eereco[i].nni > N)
                    update_nn_no_cross!(eereco, i, N, algorithm, dij_factor)
                    # update_nn_no_cross!(i, N, clusterseq, clusterseq_index, nndist,
                    #                     nndij, nni)
                end
            end
        end

        # Finally, we need to update the nearest neighbours for the new jet, checking both ways
        # (But only if there was a new jet!)
        if ijetA != ijetB
            update_nn_cross!(eereco, i, N, algorithm, dij_factor)
        end

        # ee_check_consistency(clusterseq, clusterseq_index, N, nndist, nndij, nni,
        #                      "iteration $iter")
    end

    # Return the final cluster sequence structure
    clusterseq
end
