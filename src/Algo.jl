"""
This module defines the anti-kt algorithm and similar jet reconstruction algorithms.
"""
module Algo

import ..JetReconstruction

export anti_kt, anti_kt_alt #, sequential_jet_reconstruct

function sequential_jet_reconstruct(_objects::AbstractArray{T}; p=-1, R=1, recombine=+) where T
    jets = T[] # result
    sequences = Vector{Int}[] # recombination sequences, WARNING: first index in the sequence is not necessarily the seed

    _R2 = R*R

    objects = T[]
    cyl = [[JetReconstruction.pt(_objects[1]), JetReconstruction.eta(_objects[1]), JetReconstruction.phi(_objects[1])]]; popfirst!(cyl) # cylindrical objects
    tmp_sequences = Vector{Int}[] # temporary sequences indexed according to objects
    dB = Float64[] # caches d_{iB} distances
    dp = Vector{Float64}[] # caches d_{ij} distances

    # d_{iB} distance
    function dist(i)
        if p < 0
            mp = -2p
            return @fastmath _R2/(cyl[i][1]^(mp))
        end
        @fastmath cyl[i][1]^(2p)*_R2
    end

    # d_{ij} distance
    function dist(i, j)
        cyli = cyl[i]
        cylj = cyl[j]
        d2 = cyli[2] - cylj[2] #eta
        d3 = abs(cyli[3] - cylj[3]) #phi
        if d3 > π
            d3 = 2π - d3
        end
        Δ2 = muladd(d2, d2, d3*d3)
        if p < 0 # this weird branchy "if" gives a tremendous speedup to the classic anti-kt
            mp = -2p
            if mp == 2
                m = @fastmath max(cyli[1]^(2), cylj[1]^(2))
            else
                m = @fastmath max(cyli[1]^(mp), cylj[1]^(mp))
            end
            m = @fastmath 1/m
        else
            p2 = 2p
            m = @fastmath min(cyli[1]^(p2), cylj[1]^(p2))
        end
        @fastmath m*Δ2
    end

    # a function to remove objects from the d_{ij} distance distance cache
    function remove_dp(i)
        #print(dp)
        deleteat!(dp, i)
        @inbounds for j in 1:(i-1)
            deleteat!(dp[j], i-j)
        end
        nothing
    end

    # prepare the caches
    @inbounds for i in 1:length(_objects)
        obj = _objects[i]
        push!(objects, obj)
        push!(cyl, [JetReconstruction.pt(obj), JetReconstruction.eta(obj), JetReconstruction.phi(obj)])
        push!(tmp_sequences, [i])
        push!(dB, dist(i))
    end
    @inbounds for i in 1:length(objects)
        di = Float64[]
        @inbounds for j in (i+1):length(objects)
            push!(di, dist(i, j))
        end
        push!(dp, di)
    end

    # main iteration
    while !isempty(objects)
        mindist_idx::Vector{Int64} = Int64[0, 0] # either [i, j] or [i, 0] depending on the type of the minimal found distance
        mindist = Inf
        d = 0
        @inbounds for i in 1:length(objects)
            #println("i=$i; actual i=$(tmp_sequences[i])")
            #println("\tdiB=$(dB[i]); kt2=$(cyl[i][1]^2)")
            #println("\tobj: [" * prod(["$el, " for el in _objects[i]]) *"]")
            @inbounds for j in (i+1):length(objects)
                #println("\tj=$j; actual j=$(tmp_sequences[j])")
                d = dp[i][j-i]
                #println("\t\tdp[$i][$(j-i)]=$(dp[i][j-i])")
                if d < mindist
                    mindist = d
                    mindist_idx[1] = i
                    mindist_idx[2] = j
                end
            end

            d = dB[i]
            if d < mindist
                mindist = d
                mindist_idx[1], mindist_idx[2] = i, 0
            end
        end
        #println(mindist_idx)
        #println()

        if mindist_idx[2] == 0 #if min is d_{iB}
            push!(jets, objects[mindist_idx[1]])
            push!(sequences, tmp_sequences[mindist_idx[1]])
            deleteat!(objects, mindist_idx[1])
            deleteat!(cyl, mindist_idx[1])
            deleteat!(tmp_sequences, mindist_idx[1])
            deleteat!(dB, mindist_idx[1])
            remove_dp(mindist_idx[1])
        else #if min is d_{ij}
            pseudojet = recombine(objects[mindist_idx[1]], objects[mindist_idx[2]])
            newseq = cat(tmp_sequences[mindist_idx[2]], tmp_sequences[mindist_idx[1]], dims=1) # WARNING: first index in the sequence is not necessarily the seed
            push!(objects, pseudojet)
            push!(cyl, [JetReconstruction.pt(pseudojet), JetReconstruction.eta(pseudojet), JetReconstruction.phi(pseudojet)])
            push!(tmp_sequences, newseq)
            deleteat!(objects, mindist_idx)
            deleteat!(cyl, mindist_idx)
            deleteat!(tmp_sequences, mindist_idx)
            deleteat!(dB, mindist_idx)
            remove_dp(mindist_idx[2]); remove_dp(mindist_idx[1])
            push!(dB, dist(length(objects)))
            push!(dp, Float64[])
            @inbounds for i in 1:(length(dp)-1)
                push!(dp[i], dist(i, length(dp)))
            end
        end
    end

    # return the result
    jets, sequences
end

"""
`anti_kt(objects; R=1, recombine=(x, y)->(x + y)) -> Vector, Vector{Vector{Int}}`

Runs the anti-kt jet reconstruction algorithm. `objects` can be any collection of *unique* elements.

Returns:
    `jets` - a vector of jets. Each jet is of the same type as elements in `objects`.
    `sequences` - a vector of vectors of indeces in `objects`. For all `i`, `sequences[i]` gives a sequence of indeces of objects that have been combined into the i-th jet (`jets[i]`).
"""
function anti_kt(objects; p=-1, R=1, recombine=+)
    sequential_jet_reconstruct(objects, p=p, R=R, recombine=recombine)
end

## EXPERIMENTAL ZONE
# typically sequential_jet_reconstruct_alt is used when developing an alternative (possibly better) way of reclustureing to keep the previous working version intact

function sequential_jet_reconstruct_alt(objects::AbstractArray{T}; p=-1, R=1, recombine=+) where T
    # returned values
    jets = T[] # result
    sequences = Vector{Int}[] # recombination sequences, WARNING: first index in the sequence is not necessarily the seed

    # caches
    tmp_sequences = Vector{Int}[] # temporary sequences indexed according to objects
    dB = Float64[] # caches d_{iB} distances
    dp = Vector{Float64}[] # caches d_{ij} distances

    # params
    _R2 = R*R
    _p = (round(p) == p) ? Int(p) : p # integer p if possible
    ap = abs(_p); # absolute p
    N = length(objects)

    # data
    _objects = copy(objects)
    _kt2 = JetReconstruction.pt.(_objects) .^ 2
    _phi = JetReconstruction.phi.(_objects)
    _eta = JetReconstruction.eta.(_objects)

    # d_{iB} distance times R^2
    function dist(i)
        if _p < 0
            return @fastmath _R2/(_kt2[i]^ap)
        end
        @fastmath (_kt2[i]^ap)*_R2
    end

    # d_{ij} distance times R^2
    function dist(i, j)
        deta = _eta[i] - _eta[j]
        dphi = abs(_phi[i] - _phi[j])
        if dphi > π
            dphi = 2π - dphi
        end
        Δ2 = muladd(deta, deta, dphi*dphi)
        if p < 0
            m = @fastmath 1/max(_kt2[i]^ap, _kt2[j]^ap)
        else
            m = @fastmath min(_kt2[i]^ap, _kt2[j]^ap)
        end
        @fastmath m*Δ2
    end

    # function to swap two jets in all the arrays tmp_sequences, _objects, _kt2, _phi, _eta, dB
    function swap_pseudojets(a, b)
        a, b = @fastmath minmax(a, b)

        tmp_sequences[a], tmp_sequences[b] = tmp_sequences[b], tmp_sequences[a]
        _objects[a], _objects[b] = _objects[b], _objects[a]
        _phi[a], _phi[b] = _phi[b], _phi[a]
        _eta[a], _eta[b] = _eta[b], _eta[a]
        _kt2[a], _kt2[b] = _kt2[b], _kt2[a]
        dB[a], dB[b] = dB[b], dB[a]

        # swap dp
        @inbounds for k in 1:(a-1)
            jmk = b-k; imk = a-k
            dp[k][jmk], dp[k][imk] = dp[k][imk], dp[k][jmk] # dist(k, b), dist(k, a)
        end
        @inbounds for k in (a+1):(b-1)
            jmk = b-k; kmi = k-a
            dp[a][kmi], dp[k][jmk] = dp[k][jmk], dp[a][kmi] # dist(k, a), dist(k, b)
        end
        @inbounds for k in (b+1):N
            kmj = k-b; kmi = k-a
            dp[b][kmj], dp[a][kmi] = dp[a][kmj], dp[b][kmi] # dist(k, b), dist(k, a)
        end
    end

    # prepare the caches
    @inbounds for i in 1:N
        obj = _objects[i]
        push!(tmp_sequences, [i])
        push!(dB, dist(i))

        di = Float64[]
        @inbounds for j in (i+1):N
            push!(di, dist(i, j))
        end
        push!(dp, di)
    end

    # main iteration
    while N != 0
        mindist_idx::Vector{Int64} = Int64[0, 0] # either [i, j] or [i, 0] depending on the type of the minimal found distance
        mindist = Inf
        d = 0
        @inbounds for i in 1:N
            #println("i=$i; actual i=$(tmp_sequences[i])")
            #println("\tdiB=$(dB[i]); kt2=$(_kt2[i])")
            #println("\tobj: [" * prod(["$el, " for el in _objects[i]]) *"]")
            @inbounds for j in (i+1):N
                #println("\tj=$j; actual j=$(tmp_sequences[j])")
                d = dp[i][j-i]
                #print("\t\tdp[$i][$(j-i)]=$(dp[i][j-i])")
                #(d != dp[i][j-i]) && print("; actual d=$d")
                #println()
                if d < mindist
                    mindist = d
                    mindist_idx[1], mindist_idx[2] = i, j # j > i
                end
            end

            d = dB[i]
            #(d != dB[i]) && println("\t\tdB[$i]=$(dB[i]), actual d=$d")
            if d < mindist
                mindist = d
                mindist_idx[1], mindist_idx[2] = i, 0
            end
        end
        #println(mindist_idx)
        #println()

        if mindist_idx[2] == 0 #if min is d_{iB}
            # save the jet and the corresponding recombination sequence, remove it from the current list of objects
            push!(jets, _objects[mindist_idx[1]])
            push!(sequences, tmp_sequences[mindist_idx[1]])

            swap_pseudojets(mindist_idx[1], N) # swap this jet with the last one, which gets dropped
        else #if min is d_{ij}
            # put the new combined pseudojet into the index i (mindist_idx[1]), remove the second point at j (mindist_idx[2]), recompute the distances where needed as well as kt2, phi, and eta
            _objects[mindist_idx[1]] = recombine(_objects[mindist_idx[1]], _objects[mindist_idx[2]])

            for x in tmp_sequences[mindist_idx[2]] # WARNING: first index in the sequence is not necessarily the seed
                push!(tmp_sequences[mindist_idx[1]], x)
            end
            tmp_sequences[mindist_idx[2]] = Int[]

            swap_pseudojets(mindist_idx[2], N) # swap this jet with the last one, which gets dropped

            _kt2[mindist_idx[1]] = JetReconstruction.pt(_objects[mindist_idx[1]])^2
            _phi[mindist_idx[1]] = JetReconstruction.phi(_objects[mindist_idx[1]])
            _eta[mindist_idx[1]] = JetReconstruction.eta(_objects[mindist_idx[1]])
            # update dB after we have the kt2, phi and eta values
            dB[mindist_idx[1]] = dist(mindist_idx[1])
            # update dp
            @inbounds for k in 1:(mindist_idx[1]-1)
                dp[k][mindist_idx[1]-k] = dist(k, mindist_idx[1])
            end
            @inbounds for k in (mindist_idx[1]+1):N
                dp[mindist_idx[1]][k-mindist_idx[1]] = dist(k, mindist_idx[1])
            end
        end
        N -= 1 # lazy drop the last element
    end

    # return the result
    jets, sequences
end

function anti_kt_alt(objects; p=-1, R=1, recombine=+)
    sequential_jet_reconstruct_alt(objects, p=p, R=R, recombine=recombine)
end

end
