"""
This module defines the anti-kt algorithm and similar jet reconstruction algorithms.
"""
module Algo

import ..JetReconstruction

# TODO: recombination sequences (fix the seeds?), choose better data structures, add better documentation
# TODO: function reversed_kt(objects; R=1) end

export anti_kt, anti_kt_alt #, sequential_jet_reconstruct

function sequential_jet_reconstruct(_objects::AbstractArray{T}; p=-1, R=1, recombine=+) where T

    jets = T[] # result
    sequences = Vector{Int}[] # recombination sequences, WARNING: first index in the sequence is not necessarily the seed

    objects = T[]
    cyl = [[JetReconstruction.pt(_objects[1]), JetReconstruction.eta(_objects[1]), JetReconstruction.phi(_objects[1])]]; popfirst!(cyl) # cylindrical objects
    tmp_sequences = Vector{Int}[] # temporary sequences indexed according to objects
    dB = Float64[] # caches d_{iB} distances
    dp = Vector{Float64}[] # caches d_{ij} distances

    # d_{iB} distance
    function dist(i)
        if p < 0
            mp = -2p
            return @fastmath 1/(cyl[i][1]^(mp))
        end
        @fastmath cyl[i][1]^(2p)
    end

    # d_{ij} distance
    function dist(i, j)
        cyli = cyl[i]
        cylj = cyl[j]
        d2 = cyli[2] - cylj[2]
        d3 = cyli[3] - cylj[3]
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
        @fastmath m*Δ2/(R*R)
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
            @inbounds for j in (i+1):length(objects)
                d = dp[i][j-i]
                if d < mindist
                    mindist = d
                    mindist_idx[1], mindist_idx[2] = i, j
                end
            end

            d = dB[i]
            if d < mindist
                mindist = d
                mindist_idx[1], mindist_idx[2] = i, 0
            end
        end

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
# alternative funcitons defined to try other approaches
# and test them together with the current one

function sequential_jet_reconstruct_alt(_objects::AbstractArray{T}; p=-1, R=1, recombine=+) where T
    jets = T[] # result
    sequences = Vector{Int}[] # recombination sequences, WARNING: first index in the sequence is not necessarily the seed

    objects = T[]
    cyl = [[JetReconstruction.pt(_objects[1]), JetReconstruction.eta(_objects[1]), JetReconstruction.phi(_objects[1])]]; popfirst!(cyl) # cylindrical objects
    tmp_sequences = Vector{Int}[] # temporary sequences indexed according to objects

    @inbounds for i in 1:length(_objects)
        obj = _objects[i]
        push!(objects, obj)
        push!(cyl, [JetReconstruction.pt(obj), JetReconstruction.eta(obj), JetReconstruction.phi(obj)])
        push!(tmp_sequences, [i])
    end

    # d_{ij}
    function dist(i, j)
        cyli = cyl[i]
        cylj = cyl[j]
        d2 = cyli[2] - cylj[2]
        d3 = cyli[3] - cylj[3]
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
        @fastmath m*Δ2/(R*R)
    end

    # d_{iB}
    function dist(i)
        if p < 0
            mp = -2p
            return @fastmath 1/(cyl[i][1]^(mp))
        end
        @fastmath cyl[i][1]^(2p)
    end

    while !isempty(objects)
        mindist_idx::Vector{Int64} = Int64[0, 0] # either [i, j] or [i, 0] depending on the type of the minimal found distance
        mindist = Inf
        d = 0
        @inbounds for i in 1:length(objects)
            #println("i=$i; original: i$(tmp_sequences[i])")
            @inbounds for j in (i+1):length(objects)
                #println("\t j=$j; original: j$(tmp_sequences[j])")
                d = dist(i, j)
                #println("\t d = dist($i,$j) = $d")
                if d < mindist
                    #println("\t $d < $mindist")
                    mindist = d
                    mindist_idx[1], mindist_idx[2] = i, j
                end
            end

            d = dist(i)
            #println("d = dist($i) = $d")
            if d < mindist
                #println("$d < $mindist")
                mindist = d
                mindist_idx[1], mindist_idx[2] = i, 0
            end
        end
        #println()
        if mindist_idx[2] == 0 #if min is d_{iB}
            push!(jets, objects[mindist_idx[1]])
            push!(sequences, tmp_sequences[mindist_idx[1]])
            deleteat!(objects, mindist_idx[1])
            deleteat!(cyl, mindist_idx[1])
            deleteat!(tmp_sequences, mindist_idx[1])
        else #if min is d_{ij}
            pseudojet = recombine(objects[mindist_idx[1]], objects[mindist_idx[2]])
            newseq = cat(tmp_sequences[mindist_idx[2]], tmp_sequences[mindist_idx[1]], dims=1) # WARNING: first index in the sequence is not necessarily the seed
            push!(objects, pseudojet)
            push!(cyl, [JetReconstruction.pt(pseudojet), JetReconstruction.eta(pseudojet), JetReconstruction.phi(pseudojet)])
            push!(tmp_sequences, newseq)
            deleteat!(objects, mindist_idx)
            deleteat!(cyl, mindist_idx)
            deleteat!(tmp_sequences, mindist_idx)
        end
    end

    jets, sequences
end

function anti_kt_alt(objects; p=-1, R=1, recombine=+)
    sequential_jet_reconstruct_alt(objects, p=p, R=R, recombine=recombine)
end

end
