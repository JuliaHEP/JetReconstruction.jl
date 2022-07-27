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
        Δ = muladd(d2, d2, d3*d3)
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
        @fastmath m*Δ/(R*R)
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
        mindist_idx::Vector{Int64} = Int64[1] # either [j, i] or [i] depending on the type of the minimal found distance
        mindist = Inf
        @inbounds for i in 1:length(objects)
            d = dist(i)
            if d <= mindist
                mindist = d
                mindist_idx = Int64[i]
            end
            @inbounds for j in 1:(i-1)
                d = dist(i, j)
                if d <= mindist
                    mindist = d
                    mindist_idx = Int64[j, i]
                end
            end
        end

        if length(mindist_idx) == 1 #if min is d_{iB}
            push!(jets, objects[mindist_idx[1]])
            push!(sequences, tmp_sequences[mindist_idx[1]])
        else #if min is d_{ij}
            pseudojet = recombine(objects[mindist_idx[1]], objects[mindist_idx[2]])
            newseq = cat(tmp_sequences[mindist_idx[2]], tmp_sequences[mindist_idx[1]], dims=1) # WARNING: first index in the sequence is not necessarily the seed
            push!(objects, pseudojet)
            push!(cyl, [JetReconstruction.pt(pseudojet), JetReconstruction.eta(pseudojet), JetReconstruction.phi(pseudojet)])
            push!(tmp_sequences, newseq)
        end
        deleteat!(objects, mindist_idx)
        deleteat!(cyl, mindist_idx)
        deleteat!(tmp_sequences, mindist_idx)
    end

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

# d_{ij}
function dist_alt(i, j, cyl, p, R)
    Δ = (cyl[i].eta - cyl[j].eta)^2 + (cyl[i].phi - cyl[j].phi)^2
    min(cyl[i].pt^(2p), cyl[j].pt^(2p))*Δ/(R^2)
end

# d_{iB}
function dist_alt(i, cyl, p, R)
    cyl[i].pt^(2p)
end

struct JetData{K, Eta, Phi}
    pt::K
    eta::Eta
    phi::Phi
    seq::Vector{Int}
end

function sequential_jet_reconstruct_alt!(objects::AbstractArray{T}; p=-1, R=1, recombine=+) where T
    jets = T[] # result
    sequences = Vector{Int}[] # recombination sequences, WARNING: first index in the sequence is not necessarily the seed
    pseudojets = [JetData(JetReconstruction.pt(objects[i]), JetReconstruction.eta(objects[i]), JetReconstruction.phi(objects[i]), [i]) for i in 1:length(objects)]

    while !isempty(objects)
        mindist_idx::Vector{Int64} = Int64[1, 0] # either [j, i] or [i, 0] depending on the type of the minimal found distance
        mindist = Inf
        for i in 1:length(objects)
            d = dist_alt(i, pseudojets, p, R)
            if d <= mindist
                mindist = d
                mindist_idx[1] = i
                mindist_idx[2] = 0
            end
            for j in 1:(i-1)
                d = dist_alt(i, j, pseudojets, p, R)
                if d <= mindist
                    mindist = d
                    mindist_idx[1] = j
                    mindist_idx[2] = i
                end
            end
        end

        if mindist_idx[2] == 0 #if min is d_{iB}
            push!(jets, objects[mindist_idx[1]])
            push!(sequences, pseudojets[mindist_idx[1]].seq)
            deleteat!(objects, mindist_idx[1])
            deleteat!(pseudojets, mindist_idx[1])
        else #if min is d_{ij}
            pseudojet = recombine(objects[mindist_idx[1]], objects[mindist_idx[2]])
            newseq = cat(pseudojets[mindist_idx[2]].seq, pseudojets[mindist_idx[1]].seq, dims=1) # WARNING: first index in the sequence is not necessarily the seed
            push!(objects, pseudojet)
            push!(pseudojets, JetData(JetReconstruction.pt(pseudojet), JetReconstruction.eta(pseudojet), JetReconstruction.phi(pseudojet), newseq))
            deleteat!(objects, mindist_idx)
            deleteat!(pseudojets, mindist_idx)
        end
    end

    jets, sequences
end

function anti_kt_alt(objects; p=-1, R=1, recombine=+)
    new_objects = [obj for obj in objects] # copies & converts to Vector
    sequential_jet_reconstruct_alt!(new_objects, p=p, R=R, recombine=recombine)
end

#=
function sequential_jet_reconstruct!(objects::AbstractArray{T}; p=-1, R=1, recombine=+) where T
    jets = T[] # result
    sequences = Vector{Int}[] # recombination sequences, WARNING: first index in the sequence is not necessarily the seed
    cyl = [[JetReconstruction.pt(obj), JetReconstruction.eta(obj), JetReconstruction.phi(obj)] for obj in objects] # cylindrical objects
    tmp_sequences = Vector{Int}[[i] for i in 1:length(objects)] # temporary sequences indexed according to objects

    # d_{ij}
    function dist(i, j)
        cyli = cyl[i]
        cylj = cyl[j]
        d2 = cyli[2] - cylj[2]
        d3 = cyli[3] - cylj[3]
        Δ = d2*d2 + d3*d3
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
        @fastmath m*Δ/(R*R)
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
        mindist_idx::Vector{Int64} = Int64[1] # either [j, i] or [i] depending on the type of the minimal found distance
        mindist = Inf
        @inbounds for i in 1:length(objects)
            d = dist(i)
            if d <= mindist
                mindist = d
                mindist_idx = Int64[i]
            end
            @inbounds for j in 1:(i-1)
                d = dist(i, j)
                if d <= mindist
                    mindist = d
                    mindist_idx = Int64[j, i]
                end
            end
        end

        if length(mindist_idx) == 1 #if min is d_{iB}
            push!(jets, objects[mindist_idx[1]])
            push!(sequences, tmp_sequences[mindist_idx[1]])
        else #if min is d_{ij}
            pseudojet = recombine(objects[mindist_idx[1]], objects[mindist_idx[2]])
            newseq = cat(tmp_sequences[mindist_idx[2]], tmp_sequences[mindist_idx[1]], dims=1) # WARNING: first index in the sequence is not necessarily the seed
            push!(objects, pseudojet)
            push!(cyl, [JetReconstruction.pt(pseudojet), JetReconstruction.eta(pseudojet), JetReconstruction.phi(pseudojet)])
            push!(tmp_sequences, newseq)
        end
        deleteat!(objects, mindist_idx)
        deleteat!(cyl, mindist_idx)
        deleteat!(tmp_sequences, mindist_idx)
    end

    jets, sequences
end
=#

end
