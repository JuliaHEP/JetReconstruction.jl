"""
This module defines the anti-kt algorithm and similar jet reconstruction algorithms.
"""
module Algo

import ..JetReconstruction

export anti_kt, anti_kt_alt #, sequential_jet_reconstruct

function sequential_jet_reconstruct(objects::AbstractArray{T}; p=-1, R=1, recombine=+) where T
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
        deta = abs(_eta[i] - _eta[j])
        dphi = abs(_phi[i] - _phi[j])
        if dphi > π
            dphi = 2π - dphi
        end
        Δ2 = muladd(deta, deta, dphi*dphi)
        if p < 0
            return @fastmath Δ2/max(_kt2[i]^ap, _kt2[j]^ap)
        end
        @fastmath Δ2*min(_kt2[i]^ap, _kt2[j]^ap)
    end

    # a function to remove objects from the d_{ij} distance cache
    function remove_dp(i)
        #print(dp)
        deleteat!(dp, i)
        @inbounds for j in 1:(i-1)
            deleteat!(dp[j], i-j)
        end
        nothing
    end

    # prepare the caches
    @inbounds for i in 1:N
        obj = _objects[i]
        push!(tmp_sequences, [i])
        push!(dB, dist(i))

        di = fill(0.0, N-i)
        @inbounds for j in (i+1):N
            di[j-i] = dist(i, j)
        end
        push!(dp, di)
    end

    # main iteration
    while N != 0
        mindist_idx::Vector{Int64} = Int64[0, 0] # either [i, j] or [i, 0] depending on the type of the minimal found distance
        mindist::Float64 = Inf
        d::Float64 = 0.0
        @inbounds for i in 1:N
            @inbounds for j in 1:(N-i)
                d = dp[i][j] # dist(i, j)
                if d < mindist
                    mindist = d
                    # j > i
                    mindist_idx[1] = i
                    mindist_idx[2] = j+i
                end
            end

            d = dB[i] # dist(i)
            if d < mindist
                mindist = d
                mindist_idx[1] = i
                mindist_idx[2] = 0
            end
        end

        if mindist_idx[2] == 0 #if min is d_{iB}
            push!(jets, _objects[mindist_idx[1]])
            push!(sequences, tmp_sequences[mindist_idx[1]])
            deleteat!(_objects, mindist_idx[1])
            deleteat!(_phi, mindist_idx[1])
            deleteat!(_eta, mindist_idx[1])
            deleteat!(_kt2, mindist_idx[1])
            deleteat!(tmp_sequences, mindist_idx[1])
            deleteat!(dB, mindist_idx[1])
            remove_dp(mindist_idx[1])
            N -= 1
        else #if min is d_{ij}
            # put the new recombined object to index mindist_idx[1] and remove mindist_idx[2];
            _objects[mindist_idx[1]] = recombine(_objects[mindist_idx[1]], _objects[mindist_idx[2]])

            for x in tmp_sequences[mindist_idx[2]] # WARNING: first index in the sequence is not necessarily the seed
                push!(tmp_sequences[mindist_idx[1]], x)
            end
            _phi[mindist_idx[1]] = JetReconstruction.phi(_objects[mindist_idx[1]])
            _eta[mindist_idx[1]] = JetReconstruction.eta(_objects[mindist_idx[1]])
            _kt2[mindist_idx[1]] = JetReconstruction.pt(_objects[mindist_idx[1]])^2

            deleteat!(_objects, mindist_idx[2])
            deleteat!(_phi, mindist_idx[2])
            deleteat!(_eta, mindist_idx[2])
            deleteat!(_kt2, mindist_idx[2])
            deleteat!(tmp_sequences, mindist_idx[2])
            deleteat!(dB, mindist_idx[2])
            remove_dp(mindist_idx[2])

            N -= 1

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
    dp = Vector{Float64}[] # caches d_{ij} distances between points
    dm = Vector{Bool}[] # mask for the d_{ij} distances between points

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
        deta = abs(_eta[i] - _eta[j])
        dphi = abs(_phi[i] - _phi[j])
        if dphi > π
            dphi = 2π - dphi
        end
        Δ2 = @fastmath muladd(deta, deta, dphi*dphi)
        if Δ2 >= _R2 # exit earlier if the distance is irrelevant
            dm[i][j-i] = false
            return Inf
        end
        if p < 0
            return @fastmath Δ2/max(_kt2[i]^ap, _kt2[j]^ap)
        end
        @fastmath Δ2*min(_kt2[i]^ap, _kt2[j]^ap)
    end

    # a function to remove objects from dp and dm arrays
    function remove_from_cache(d, i)
        deleteat!(d, i)
        @inbounds for j in 1:(i-1)
            deleteat!(d[j], i-j)
        end
        nothing
    end

    # prepare the caches
    @inbounds for i in 1:N
        obj = _objects[i]
        push!(tmp_sequences, [i])
        push!(dB, dist(i))

        dmi = fill(true, N-i)
        push!(dm, dmi)

        di = fill(0.0, N-i)
        @inbounds for j in (i+1):N
            di[j-i] = dist(i, j)
        end
        push!(dp, di)
    end

    # main iteration
    while N != 0
        mindist_idx::Vector{Int64} = Int64[0, 0] # either [i, j] or [i, 0] depending on the type of the minimal found distance
        mindist::Float64 = Inf
        d::Float64 = 0.0
        @inbounds for i in 1:N
            @inbounds for j in 1:(N-i)
                if !dm[i][j]
                    continue
                end
                d = dp[i][j] # dist(i, j)
                if d < mindist
                    mindist = d
                    # j > i
                    mindist_idx[1] = i
                    mindist_idx[2] = j+i
                end
            end

            d = dB[i] # dist(i)
            if d < mindist
                mindist = d
                mindist_idx[1] = i
                mindist_idx[2] = 0
            end
        end

        if mindist_idx[2] == 0 #if min is d_{iB}
            push!(jets, _objects[mindist_idx[1]])
            push!(sequences, tmp_sequences[mindist_idx[1]])
            deleteat!(_objects, mindist_idx[1])
            deleteat!(_phi, mindist_idx[1])
            deleteat!(_eta, mindist_idx[1])
            deleteat!(_kt2, mindist_idx[1])
            deleteat!(tmp_sequences, mindist_idx[1])
            deleteat!(dB, mindist_idx[1])
            remove_from_cache(dp, mindist_idx[1])
            remove_from_cache(dm, mindist_idx[1])
            N -= 1
        else #if min is d_{ij}
            # put the new recombined object to index mindist_idx[1] and remove mindist_idx[2];
            _objects[mindist_idx[1]] = recombine(_objects[mindist_idx[1]], _objects[mindist_idx[2]])

            for x in tmp_sequences[mindist_idx[2]] # WARNING: first index in the sequence is not necessarily the seed
                push!(tmp_sequences[mindist_idx[1]], x)
            end
            _phi[mindist_idx[1]] = JetReconstruction.phi(_objects[mindist_idx[1]])
            _eta[mindist_idx[1]] = JetReconstruction.eta(_objects[mindist_idx[1]])
            _kt2[mindist_idx[1]] = JetReconstruction.pt(_objects[mindist_idx[1]])^2

            deleteat!(_objects, mindist_idx[2])
            deleteat!(_phi, mindist_idx[2])
            deleteat!(_eta, mindist_idx[2])
            deleteat!(_kt2, mindist_idx[2])
            deleteat!(tmp_sequences, mindist_idx[2])
            deleteat!(dB, mindist_idx[2])
            remove_from_cache(dp, mindist_idx[2])
            remove_from_cache(dm, mindist_idx[2])

            N -= 1

            # update dB after we have the kt2, phi and eta values
            dB[mindist_idx[1]] = dist(mindist_idx[1])
            # update dp (dm gets updated automatically)
            @inbounds for k in 1:(mindist_idx[1]-1)
                dp[k][mindist_idx[1]-k] = dist(k, mindist_idx[1])
            end
            @inbounds for k in (mindist_idx[1]+1):N
                dp[mindist_idx[1]][k-mindist_idx[1]] = dist(mindist_idx[1], k)
            end
        end
    end

    # return the result
    jets, sequences
end

"""
Not for usage. Use `anti_kt` instead. Correctness is not guaranteed.
"""
function anti_kt_alt(objects; p=-1, R=1, recombine=+)
    sequential_jet_reconstruct_alt(objects, p=p, R=R, recombine=recombine)
end

end
