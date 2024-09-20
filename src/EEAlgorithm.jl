# Use these constants whenever we need to set a large value for the distance or
# metric distance
const large_distance = 16.0 # = 4^2
const large_dij = 1.0e6

"""
    angular_distance(eereco, i, j) -> Float64

Calculate the angular distance between two jets using the formula (1 - cos(θ)).

# Arguments
- `eereco`: The array of `EERecoJet` objects.
- `i`: The first jet.
- `j`: The second jet.

# Returns
- `Float64`: The angular distance between `i` and `j`, which is ``1 - cos\theta``.
"""
@inline function angular_distance(eereco, i, j)
    @muladd 1.0 - eereco[i].nx * eereco[j].nx - eereco[i].ny * eereco[j].ny - eereco[i].nz * eereco[j].nz
end

"""
    dij_dist(eereco, i, j, dij_factor)

Calculate the dij distance between two ``e^+e^-``jets.

# Arguments
- `eereco`: The array of `EERecoJet` objects.
- `i`: The first jet.
- `j`: The second jet.
- `dij_factor`: The scaling factor to multiply the dij distance by.

# Returns
- The dij distance between `i` and `j`.
"""
@inline function dij_dist(eereco, i, j, dij_factor)
    # Calculate the dij distance for jet i
    min(eereco[i].E2p, eereco[j].E2p) * dij_factor * eereco[i].nndist
end

function get_angular_nearest_neighbours!(eereco, algorithm, dij_factor)
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
            # @muladd this_nndist = 1.0 - eereco[i].nx * eereco[j].nx - eereco[i].ny * eereco[j].ny - eereco[i].nz * eereco[j].nz
            better_nndist_i = this_nndist < eereco[i].nndist
            eereco.nndist[i] = better_nndist_i ? this_nndist : eereco.nndist[i]
            eereco.nni[i] = better_nndist_i ? j : eereco.nni[i]
            better_nndist_j = this_nndist < eereco[j].nndist
            eereco.nndist[j] = better_nndist_j ? this_nndist : eereco.nndist[i]
            eereco.nni[j] = better_nndist_j ? i : eereco.nni[j]
            # if this_nndist < eereco[i].nndist
            #     eereco.nndist[i] = this_nndist
            #     eereco.nni[i] = j
            # end
            # if this_nndist < eereco[j].nndist
            #     eereco.nndist[j] = this_nndist
            #     eereco.nni[j] = i
            # end
        end
    end
    # Nearest neighbour dij distance
    @inbounds for i in 1:N
        eereco.dijdist[i] = min(eereco[i].E2p, eereco[eereco[i].nni].E2p) * dij_factor * eereco[i].nndist
    end
    # For the EEKt algorithm, we need to check the beam distance as well
    # This is structured to only check for EEKt once!
    if algorithm == JetAlgorithm.EEKt
        @inbounds for i in 1:N
            beam_close = eereco[i].E2p < eereco[i].dijdist
            eereco.dijdist[i] = beam_close ? eereco[i].E2p : eereco.dijdist[i]
            eereco.nni[i] = beam_close ? 0 : eereco.nni[i]
            # if eereco[i].E2p < eereco[i].dijdist
            #     eereco.dijdist[i] = eereco.E2p[i]
            #     eereco.nni[i] = 0
            # end
        end
    end
end

# Update the nearest neighbour for jet i, w.r.t. all other active jets
function update_nn_no_cross!(eereco, i, N, algorithm, dij_factor)
    eereco.nndist[i] = large_distance
    eereco.nni[i] = i
    @inbounds for j in 1:N
        if j != i
            this_nndist = angular_distance(eereco, i, j)
            # @muladd this_nndist = 1.0 - eereco[i].nx * eereco[j].nx - eereco[i].ny * eereco[j].ny - eereco[i].nz * eereco[j].nz
            better_nndist_i = this_nndist < eereco[i].nndist
            eereco.nndist[i] = better_nndist_i ? this_nndist : eereco.nndist[i]
            eereco.nni[i] = better_nndist_i ? j : eereco.nni[i]
            # if this_nndist < eereco[i].nndist
            #     eereco.nndist[i] = this_nndist
            #     eereco.nni[i] = j
            # end
        end
    end
    eereco.dijdist[i] = min(eereco[i].E2p, eereco[eereco[i].nni].E2p) * dij_factor * eereco[i].nndist
    if algorithm == JetAlgorithm.EEKt
        beam_close = eereco[i].E2p < eereco[i].dijdist
        eereco.dijdist[i] = beam_close ? eereco[i].E2p : eereco.dijdist[i]
        eereco.nni[i] = beam_close ? 0 : eereco.nni[i]
        # if eereco[i].E2p < eereco[i].dijdist
        #     eereco.dijdist[i] = eereco[i].E2p
        #     eereco.nni[i] = 0
        # end
    end
end

function update_nn_cross!(eereco, i, N, algorithm, dij_factor)
    # Update the nearest neighbour for jet i, w.r.t. all other active jets
    # also doing the cross check for the other jet
    eereco.nndist[i] = large_distance
    eereco.nni[i] = i
    @inbounds for j in 1:N
        if j != i
            this_nndist = angular_distance(eereco, i, j)
            # @muladd this_nndist = 1.0 - eereco[i].nx * eereco[j].nx - eereco[i].ny * eereco[j].ny - eereco[i].nz * eereco[j].nz
            better_nndist_i = this_nndist < eereco[i].nndist
            eereco.nndist[i] = better_nndist_i ? this_nndist : eereco.nndist[i]
            eereco.nni[i] = better_nndist_i ? j : eereco.nni[i]
            # if this_nndist < eereco[i].nndist
            #     eereco.nndist[i] = this_nndist
            #     eereco.nni[i] = j
            # end
            if this_nndist < eereco[j].nndist
                eereco.nndist[j] = this_nndist
                eereco.nni[j] = i
                # j will not be revisited, so update metric distance here
                eereco.dijdist[j] = min(eereco[i].E2p, eereco[j].E2p) * dij_factor * eereco[j].nndist
                if algorithm == JetAlgorithm.EEKt
                    if eereco[j].E2p < eereco[j].dijdist
                        eereco.dijdist[j] = eereco[j].E2p
                        eereco.nni[j] = 0
                    end
                end
            end
        end
    end
    eereco.dijdist[i] = min(eereco[i].E2p, eereco[eereco[i].nni].E2p) * dij_factor * eereco[i].nndist
    if algorithm == JetAlgorithm.EEKt
        beam_close = eereco[i].E2p < eereco[i].dijdist
        eereco.dijdist[i] = beam_close ? eereco[i].E2p : eereco.dijdist[i]
        eereco.nni[i] = beam_close ? 0 : eereco.nni[i]
        # if eereco[i].E2p < eereco[i].dijdist
        #     eereco.dijdist[i] = eereco[i].E2p
        #     eereco.nni[i] = 0
        # end
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

    # return clusterseq

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

        # println("N $N, Iteration $iter: dij_min = $dij_min, ijetA = $ijetA, ijetB = $ijetB; njets = $(length(clusterseq.jets))")

        # Normalise the dij_min
        # dij_min *= dij_factor

        # Now we check if there is a "beam" merge possibility
        if ijetB == 0
            # We have an EEKt beam merge
            ijetB = ijetA
            # println("Beam merge: dij_min = $dij_min, ijetA = $ijetA, ijetB = $ijetB")
            add_step_to_history!(clusterseq,
                                 clusterseq.jets[eereco[ijetA].index]._cluster_hist_index,
                                 BeamJet, Invalid, dij_min)
        elseif N == 1
            # We have a single jet left
            ijetB = ijetA
            # println("Final jet: dij_min = $dij_min, ijetA = $ijetA, ijetB = $ijetB")
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
            eereco.index[ijetA] = newjet_k
            eereco.nni[ijetA] = 0
            eereco.nndist[ijetA] = R2
            eereco.nx[ijetA] = nx(merged_jet)
            eereco.ny[ijetA] = ny(merged_jet)
            eereco.nz[ijetA] = nz(merged_jet)
            eereco.E2p[ijetA] = energy(merged_jet)^(2p)
        end

        # Squash step - copy the final jet's compact data into the jetB slot
        # unless we are at the end of the array, in which case do nothing
        if ijetB != N
            eereco.index[ijetB] = eereco.index[N]
            eereco.nni[ijetB] = eereco.nni[N]
            eereco.nndist[ijetB] = eereco.nndist[N]
            eereco.dijdist[ijetB] = eereco.dijdist[N]
            eereco.nx[ijetB] = eereco.nx[N]
            eereco.ny[ijetB] = eereco.ny[N]
            eereco.nz[ijetB] = eereco.nz[N]
            eereco.E2p[ijetB] = eereco.E2p[N]
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
            update_nn_cross!(eereco, ijetA, N, algorithm, dij_factor)
        end

        # ee_check_consistency(clusterseq, clusterseq_index, N, nndist, nndij, nni,
        #                      "iteration $iter")
        # ee_check_consistency(clusterseq, eereco, N)
    end

    # Return the final cluster sequence structure
    clusterseq
end
