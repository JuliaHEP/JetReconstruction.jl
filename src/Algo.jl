Base.@propagate_inbounds function _dist(i, j, _eta, _phi)
    deta = _eta[i] - _eta[j]
    dphi = abs(_phi[i] - _phi[j])
    dphi = ifelse(dphi > pi, 2pi - dphi, dphi)
    muladd(deta, deta, dphi * dphi)
end

# d_{ij} distance with i's NN (times R^2)
Base.@propagate_inbounds function _dij(i, _kt2, _nn, _nndist)
    j = _nn[i]
    d = _nndist[i]
    d * min(_kt2[i], _kt2[j])
end

# finds new nn for i and checks everyone additionally
Base.@propagate_inbounds function _upd_nn_crosscheck!(i::Int, from::Int, to::Int, _eta, _phi, _R2, _nndist, _nn)
    nndist = _R2
    nn = i
    @inbounds @simd for j in from:to
        Δ2 = _dist(i, j, _eta, _phi)
        if Δ2 < nndist
            nn = j
            nndist = Δ2
        end
        if Δ2 < _nndist[j]
            _nndist[j] = Δ2
            _nn[j] = i
        end
    end
    _nndist[i] = nndist
    _nn[i] = nn
end

# finds new nn for i
Base.@propagate_inbounds function _upd_nn_nocross!(i::Int, from::Int, to::Int, _eta, _phi, _R2, _nndist, _nn)
    nndist = _R2
    nn = i
    @inbounds @simd for j in from:(i-1)
        Δ2 = _dist(i, j, _eta, _phi)
        if Δ2 <= nndist
            nn = j
            nndist = Δ2
        end
    end
    @inbounds @simd for j in (i+1):to
        Δ2 = _dist(i, j, _eta, _phi)
        f = Δ2 <= nndist
        nn = ifelse(f, j, nn)
        nndist = ifelse(f, Δ2, nndist)
    end
    _nndist[i] = nndist
    _nn[i] = nn
end

# entire NN update step
# Base.@propagate_inbounds function _upd_nn_step!(i::Int, j::Int, k::Int, N::Int, Nn::Int, 
#     _kt2::Vector{Float64}, _eta::Vector{Float64}, _phi::Vector{Float64}, _R2::Float64, _nndist::Vector{Float64}, _nn::Vector{Int}, _nndij::Vector{Float64})
Base.@propagate_inbounds function _upd_nn_step!(i, j, k, N, Nn, _kt2, _eta, _phi, _R2, _nndist, _nn, _nndij)
    nnk = _nn[k]
    if nnk == i || nnk == j
        _upd_nn_nocross!(k, 1, N, _eta, _phi, _R2, _nndist, _nn) # update dist and nn
        _nndij[k] = _dij(k, _kt2, _nn, _nndist)
        nnk = _nn[k]
    end

    if j != i && k != i
        Δ2 = _dist(i, k, _eta, _phi)
        if Δ2 < _nndist[k]
            _nndist[k] = Δ2
            nnk = _nn[k] = i
            _nndij[k] = _dij(k, _kt2, _nn, _nndist)
        end

        cond = Δ2 < _nndist[i]
        _nndist[i], _nn[i] = ifelse(cond, (Δ2, k), (_nndist[i], _nn[i]))
    end

    nnk == Nn && (_nn[k] = j)
end


function sequential_jet_reconstruct(objects::AbstractArray{T}; p = -1.0, R = 1.0, recombine = +) where T
    # Bounds
    N::Int = length(objects)

    # Returned values
    jets = T[]
    sequences = Vector{Int}[] # recombination sequences, WARNING: first index in the sequence is not necessarily the seed

    # Parameters
    _R2 = R * R
    _p = (round(p) == p) ? Int(p) : p # integer p if possible

    # Data
    _objects = copy(objects)
    ## When working with LorentzVectorHEP we make sure these arrays are type stable
    _kt2::Vector{Float64} = (LorentzVectorHEP.pt.(_objects) .^ 2) .^ _p
    _phi::Vector{Float64} = LorentzVectorHEP.phi.(_objects)
    _eta::Vector{Float64} = LorentzVectorHEP.rap.(_objects)
    _nn = Vector(1:N) # nearest neighbours
    _nndist = fill(float(_R2), N) # distances to the nearest neighbour
    _sequences = Vector{Int}[[x] for x in 1:N]

    # initialize _nn
    @simd for i in 1:N
        _upd_nn_crosscheck!(i, 1, i - 1, _eta, _phi, _R2, _nndist, _nn)
    end

    # diJ table *_R2
    _nndij::Vector{Float64} = zeros(N)
    @inbounds @simd for i in 1:N
        _nndij[i] = _dij(i, _kt2, _nn, _nndist)
    end

    iteration::Int = 1
    while N != 0
        # findmin
        i = 1
        dij_min = _nndij[1]
        @inbounds @simd for k in 2:N
            cond = _nndij[k] < dij_min
            dij_min, i = ifelse(cond, (_nndij[k], k), (dij_min, i))
        end

        j::Int = _nn[i]

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
            _objects[i] = recombine(_objects[i], _objects[j])
            _phi[i] = LorentzVectorHEP.phi(_objects[i])
            _eta[i] = LorentzVectorHEP.rap(_objects[i])
            _kt2[i] = (LorentzVectorHEP.pt(_objects[i])^2)^_p

            _nndist[i] = _R2
            _nn[i] = i

            @inbounds for x in _sequences[j] # WARNING: first index in the sequence is not necessarily the seed
                push!(_sequences[i], x)
            end
        else # i == j
            # push!(jets, LorentzVectorCyl(LorentzVectorHEP.pt(_objects[i]), LorentzVectorHEP.eta(_objects[i]),
            #     LorentzVectorHEP.phi(_objects[i]), LorentzVectorHEP.mt(_objects[i])))
            push!(jets, _objects[i])
            push!(sequences, _sequences[i])
        end

        # copy jet N to j
        _objects[j] = _objects[N]

        _phi[j] = _phi[N]
        _eta[j] = _eta[N]
        _kt2[j] = _kt2[N]
        _nndist[j] = _nndist[N]
        _nn[j] = _nn[N]
        _nndij[j] = _nndij[N]

        _sequences[j] = _sequences[N]

        Nn::Int = N
        N -= 1
        iteration += 1

        # update nearest neighbours step
        @inbounds @simd for k in 1:N
            _upd_nn_step!(i, j, k, N, Nn, _kt2, _eta, _phi, _R2, _nndist, _nn, _nndij)
        end
        # @infiltrate

        _nndij[i] = _dij(i, _kt2, _nn, _nndist)
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
