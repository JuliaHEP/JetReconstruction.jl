import IfElse

@inline function _dist(i, j, _eta, _phi)
    deta = _eta[i] - _eta[j]
    dphi = abs(_phi[i] - _phi[j])
    dphi = IfElse.ifelse(dphi > pi, 2pi - dphi, dphi)
    muladd(deta, deta, dphi*dphi)
end

# d_{ij} distance with i's NN (times R^2)
function _dij(i, _kt2, _nn, _nndist)
    j = _nn[i]
    d = _nndist[i]
    d*min(_kt2[i], _kt2[j])
end

# finds new nn for i and checks everyone additionally
function _upd_nn_crosscheck!(i::Int, from::Int, to::Int, _eta, _phi, _R2, _nndist, _nn)
    nndist = _R2
    nn = i
    @inbounds for j in from:to
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
    nothing
end

# finds new nn for i
function _upd_nn_nocross!(i::Int, from::Int, to::Int, _eta, _phi, _R2, _nndist, _nn)
    nndist = _R2
    nn = i
    @turbo for j in from:(i-1)
        Δ2 = _dist(i, j, _eta, _phi)
        f = Δ2 <= nndist
        nn = IfElse.ifelse(f, j, nn)
        nndist = IfElse.ifelse(f, Δ2, nndist)
    end
    @turbo for j in (i+1):to
        Δ2 = _dist(i, j, _eta, _phi)
        f = Δ2 <= nndist
        nn = IfElse.ifelse(f, j, nn)
        nndist = IfElse.ifelse(f, Δ2, nndist)
    end
    _nndist[i] = nndist
    _nn[i] = nn
    nothing
end

# entire NN update step
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
        _nndist[i] = IfElse.ifelse(cond, Δ2, _nndist[i])
        _nn[i] = IfElse.ifelse(cond, k, _nn[i])
    end

    nnk == Nn && (_nn[k] = j)
end

function sequential_jet_reconstruct(objects::AbstractArray{T}; p=-1, R=1., recombine=+) where T
    # bounds
    N::Int = length(objects)

    # returned values
    jets = T[] # result
    sequences = Vector{Int}[] # recombination sequences, WARNING: first index in the sequence is not necessarily the seed

    # params
    _R2 = R*R
    _p = (round(p) == p) ? Int(p) : p # integer p if possible
    ap = abs(_p); # absolute p

    # data
    _objects = copy(objects)
    _kt2 = (JetReconstruction.pt.(_objects) .^ 2) .^ _p
    _phi = JetReconstruction.phi.(_objects)
    _eta = JetReconstruction.eta.(_objects)
    _nn = Vector(1:N) # nearest neighbours
    _nndist = fill(float(_R2), N) # distances to the nearest neighbour
    _sequences = Vector{Int}[[x] for x in 1:N]

    # initialize _nn
    for i in 1:N
        _upd_nn_crosscheck!(i, 1, i-1, _eta, _phi, _R2, _nndist, _nn)
    end

    # diJ table *_R2
    _nndij = zeros(N)
    for i in 1:N
        _nndij[i] = _dij(i, _kt2, _nn, _nndist)
    end

    while N != 0
        # findmin
        i = 1
        dij_min = _nndij[1]
        for k in 2:N
            if _nndij[k] < dij_min
                dij_min = _nndij[k]
                i = k
            end
        end

        j::Int = _nn[i]

        if i != j
            # swap if needed
            if j < i
                i, j = j, i
            end

            # update ith jet, replacing it with the new one
            _objects[i] = recombine(_objects[i], _objects[j])
            _phi[i] = JetReconstruction.phi(_objects[i])
            _eta[i] = JetReconstruction.eta(_objects[i])
            _kt2[i] = (JetReconstruction.pt(_objects[i])^2)^_p

            _nndist[i] = _R2
            _nn[i] = i

            for x in _sequences[j] # WARNING: first index in the sequence is not necessarily the seed
                push!(_sequences[i], x)
            end
        else # i == j
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

        # update nearest neighbours step
        for k::Int in 1:N
            _upd_nn_step!(i, j, k, N, Nn, _kt2, _eta, _phi, _R2, _nndist, _nn, _nndij)
        end

        _nndij[i] = _dij(i, _kt2, _nn, _nndist)
    end

    jets, sequences
end

"""
`anti_kt_algo(objects; R=1, recombine=(x, y)->(x + y)) -> Vector, Vector{Vector{Int}}`

Runs the anti-kt jet reconstruction algorithm. `objects` can be any collection of *unique* elements.

Returns:
    `jets` - a vector of jets. Each jet is of the same type as elements in `objects`.
    `sequences` - a vector of vectors of indeces in `objects`. For all `i`, `sequences[i]` gives a sequence of indeces of objects that have been combined into the i-th jet (`jets[i]`).
"""
function anti_kt_algo(objects; p=-1, R=1., recombine=+)
    sequential_jet_reconstruct(objects, p=p, R=R, recombine=recombine)
end

"""
`kt_algo(objects; R=1, recombine=(x, y)->(x + y)) -> Vector, Vector{Vector{Int}}`

Runs the kt jet reconstruction algorithm. `objects` can be any collection of *unique* elements.

Returns:
    `jets` - a vector of jets. Each jet is of the same type as elements in `objects`.
    `sequences` - a vector of vectors of indeces in `objects`. For all `i`, `sequences[i]` gives a sequence of indeces of objects that have been combined into the i-th jet (`jets[i]`).
"""
function kt_algo(objects; p=1, R=1., recombine=+)
    sequential_jet_reconstruct(objects, p=p, R=R, recombine=recombine)
end

## EXPERIMENTAL ZONE
# typically sequential_jet_reconstruct_alt is used when developing an alternative (possibly better) way of reclustureing to keep the previous working version intact

function sequential_jet_reconstruct_alt(objects::AbstractArray{T}; p=-1, R=1, recombine=+) where T
    # bounds
    N::Int = length(objects)

    # returned values
    jets = T[] # result
    sequences = Vector{Int}[] # recombination sequences, WARNING: first index in the sequence is not necessarily the seed

    # params
    _R2::Float64 = R*R
    _p = (round(p) == p) ? Int(p) : p # integer p if possible
    ap = abs(_p); # absolute p

    # data
    _objects = copy(objects)
    _kt2 = (JetReconstruction.pt.(_objects) .^ 2) .^ _p
    _phi = JetReconstruction.phi.(_objects)
    _eta = JetReconstruction.eta.(_objects)
    _nn = Vector(1:N) # nearest neighbours
    _nndist = fill(float(_R2), N) # distances to the nearest neighbour
    _sequences = Vector{Int}[[x] for x in 1:N]

    # initialize _nn
    for i::Int in 1:N
        _upd_nn_crosscheck!(i, 1, i-1, _eta, _phi, _R2, _nndist, _nn)
    end

    # diJ table *_R2
    _nndij = zeros(N)
    for i::Int in 1:N
        _nndij[i] = _dij(i, _kt2, _nn, _nndist)
    end

    while N != 0
        # findmin
        i::Int = 1
        dij_min = _nndij[1]
        for k::Int in 2:N
            if _nndij[k] < dij_min
                dij_min = _nndij[k]
                i = k
            end
        end

        j::Int = _nn[i]

        if i != j
            # swap if needed
            if j < i
                i, j = j, i
            end

            # update ith jet, replacing it with the new one
            _objects[i] = recombine(_objects[i], _objects[j])
            _phi[i] = JetReconstruction.phi(_objects[i])
            _eta[i] = JetReconstruction.eta(_objects[i])
            _kt2[i] = (JetReconstruction.pt(_objects[i])^2)^_p

            _nndist[i] = _R2
            _nn[i] = i

            for x in _sequences[j] # WARNING: first index in the sequence is not necessarily the seed
                push!(_sequences[i], x)
            end
        else # i == j
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

        # update nearest neighbours step
        for k::Int in 1:N
            _upd_nn_step!(i, j, k, N, Nn, _kt2, _eta, _phi, _R2, _nndist, _nn, _nndij)
        end

        _nndij[i] = _dij(i, _kt2, _nn, _nndist)
    end

    jets, sequences
end

"""
Not for usage. Use `anti_kt_algo` instead. Correctness is not guaranteed.
"""
function anti_kt_algo_alt(objects; p=-1, R=1, recombine=+)
    sequential_jet_reconstruct_alt(objects, p=p, R=R, recombine=recombine)
end
