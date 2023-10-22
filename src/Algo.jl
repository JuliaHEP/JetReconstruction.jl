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
of objects, which supports the methods pt2(), phi(), rapidity() for each element.
"""
function sequential_jet_reconstruct(objects::AbstractArray{T}; p = -1, R = 1.0, recombine = +) where T
    # Integer p if possible
    p = (round(p) == p) ? Int(p) : p 

    # We make sure these arrays are type stable - have seen issues where, depending on the values
    # returned by the methods, they can become unstable and performance degrades
    kt2_array::Vector{Float64} = pt2.(objects) .^ p
    phi_array::Vector{Float64} = phi.(objects)
    rapidity_array::Vector{Float64} = rapidity.(objects)

    objects_array = copy(objects)

    # Now call the actual reconstruction method, tuned for our internal EDM
    sequential_jet_reconstruct(objects_array=objects_array, kt2_array=kt2_array, phi_array=phi_array, 
        rapidity_array=rapidity_array, p=p, R=R, recombine=recombine)
end



function sequential_jet_reconstruct(;objects_array::AbstractArray{J}, kt2_array::Vector{F}, 
        phi_array::Vector{F}, rapidity_array::Vector{F}, p = -1, R = 1.0, recombine = +) where {J, F<:AbstractFloat}
    # Bounds
    N::Int = length(objects_array)

    # Returned values
    jets = J[]
    sequences = Vector{Int}[] # recombination sequences, WARNING: first index in the sequence is not necessarily the seed

    # Parameters
    R2 = R * R

    # Data
    nn = Vector(1:N) # nearest neighbours
    nndist = fill(float(R2), N) # distances to the nearest neighbour
    sequences = Vector{Int}[[x] for x in 1:N]

    # initialize _nn
    @simd for i in 1:N
        upd_nn_crosscheck!(i, 1, i - 1, rapidity_array, phi_array, R2, nndist, nn)
    end

    # diJ table *_R2
    nndij::Vector{Float64} = zeros(N)
    @inbounds @simd for i in 1:N
        nndij[i] = dij(i, kt2_array, nn, nndist)
    end

    iteration::Int = 1
    while N != 0
        # findmin
        i = 1
        dij_min = nndij[1]
        @inbounds @simd for k in 2:N
            cond = nndij[k] < dij_min
            dij_min, i = ifelse(cond, (nndij[k], k), (dij_min, i))
        end

        j::Int = nn[i]

        ## Needed for certain tricky debugging situations
        # if iteration==1
        #     debug_jets(_nn, _nndist, _nndij)
        # end

        if i != j
            # swap if needed
            if j < i
                i, j = j, i
            end

            # update ith jet, replacing it with the new one
            objects_array[i] = recombine(objects_array[i], objects_array[j])
            phi_array[i] = LorentzVectorHEP.phi(objects_array[i])
            rapidity_array[i] = LorentzVectorHEP.rapidity(objects_array[i])
            kt2_array[i] = LorentzVectorHEP.pt2(objects_array[i]) ^ p

            nndist[i] = R2
            nn[i] = i

            @inbounds for x in sequences[j] # WARNING: first index in the sequence is not necessarily the seed
                push!(sequences[i], x)
            end
        else # i == j
            push!(jets, objects_array[i])
            push!(sequences, sequences[i])
        end

        # copy jet N to j
        objects_array[j] = objects_array[N]

        phi_array[j] = phi_array[N]
        rapidity_array[j] = rapidity_array[N]
        kt2_array[j] = kt2_array[N]
        nndist[j] = nndist[N]
        nn[j] = nn[N]
        nndij[j] = nndij[N]

        sequences[j] = sequences[N]

        Nn::Int = N
        N -= 1
        iteration += 1

        # update nearest neighbours step
        @inbounds @simd for k in 1:N
            upd_nn_step!(i, j, k, N, Nn, kt2_array, rapidity_array, phi_array, R2, nndist, nn, nndij)
        end
        # @infiltrate

        nndij[i] = dij(i, kt2_array, nn, nndist)
    end

    jets, sequences
end

"""
`anti_kt_algo(objects; R=1, recombine=(x, y)->(x + y)) -> Vector, Vector{Vector{Int}}`

Runs the anti-kt jet reconstruction algorithm. `objects` can be any collection of *unique* elements.

Returns:
    `jets` - a vector of jets. Each jet is of the same type as elements in `objects`.
    `sequences` - a vector of vectors of indices in `objects`. For all `i`, `sequences[i]` gives a sequence of indices of objects that have been combined into the i-th jet (`jets[i]`).
"""
function anti_kt_algo(objects; R = 1.0, recombine = +)
    sequential_jet_reconstruct(objects, p = -1, R = R, recombine = recombine)
end

"""
`kt_algo(objects; R=1, recombine=(x, y)->(x + y)) -> Vector, Vector{Vector{Int}}`

Runs the kt jet reconstruction algorithm. `objects` can be any collection of *unique* elements.

Returns:
    `jets` - a vector of jets. Each jet is of the same type as elements in `objects`.
    `sequences` - a vector of vectors of indices in `objects`. For all `i`, `sequences[i]` gives a sequence of indices of objects that have been combined into the i-th jet (`jets[i]`).
"""
function kt_algo(objects; R = 1.0, recombine = +)
    sequential_jet_reconstruct(objects, p = 1, R = R, recombine = recombine)
end

"""
`cambridge_aachen_algo(objects; R=1, recombine=(x, y)->(x + y)) -> Vector, Vector{Vector{Int}}`

Runs the Cambridge/Aachen jet reconstruction algorithm. `objects` can be any collection of *unique* elements.

Returns:
    `jets` - a vector of jets. Each jet is of the same type as elements in `objects`.
    `sequences` - a vector of vectors of indices in `objects`. For all `i`, `sequences[i]` gives a sequence of indices of objects that have been combined into the i-th jet (`jets[i]`).
"""
function cambridge_aachen_algo(objects; R = 1.0, recombine = +)
    sequential_jet_reconstruct(objects, p = 0, R = R, recombine = recombine)
end
