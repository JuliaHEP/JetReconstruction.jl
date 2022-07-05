"""
This module defines the anti-kt algorithm and similar jet reconstruction algorithms.
"""
module Algo

# TODO: recombination tree, choose better data structures, add better documentation
# TODO: function reversed_kt(objects; R=1) end

export anti_kt!, anti_kt #, sequential_jet_reconstruct!, sequential_jet_reconstruct

function sequential_jet_reconstruct!(objects::AbstractArray{T}; p=-1, R=1, recombine=+) where T
    jets = T[] # result
    cyl = [[Main.pt(obj), Main.eta(obj), Main.phi(obj)] for obj in objects] # cylindrical objects SHOULD WE CALCULATE THEM HERE OR LATER? Maybe switch to StaticVector (or only if installed)

    # d_{ij}
    function dist(i, j)
        Δ = (cyl[i][2] - cyl[j][2])^2 + (cyl[i][3] - cyl[j][3])^2
        min(cyl[i][1]^(2p), cyl[j][1]^(2p))*Δ/(R^2)
    end

    # d_{iB}
    function dist(i)
        cyl[i][1]^(2p)
    end

    while !isempty(objects)
        mindist_idx::Vector{Int64} = Int64[1] # either [j, i] or [i] depending on the type of the minimal found distance
        mindist = Inf
        for i in 1:length(objects)
            d = dist(i)
            if d < mindist
                mindist = d
                mindist_idx = Int64[i]
            end
            for j in 1:(i-1)
                d = dist(i, j)
                if d < mindist
                    mindist = d
                    mindist_idx = Int64[j, i]
                end
            end
        end

        if length(mindist_idx) == 1 #if min is d_{iB}
            push!(jets, objects[mindist_idx[1]])
        else #if min is d_{ij}
            pseudojet = recombine(objects[mindist_idx[1]], objects[mindist_idx[2]])
            push!(objects, pseudojet)
            push!(cyl, [Main.pt(pseudojet), Main.eta(pseudojet), Main.phi(pseudojet)])
        end
        deleteat!(objects, mindist_idx)
        deleteat!(cyl, mindist_idx)
    end

    jets#, tree
end

function sequential_jet_reconstruct(objects; p=-1, R=1, recombine=+)
    new_objects = [obj for obj in objects] # copies & converts to Vector
    sequential_jet_reconstruct!(new_objects, p=p, R=R, recombine=recombine)
end

"""
`anti_kt!(objects::AbstractArray{T} where T; R=1, recombine=(x, y)->(x + y)) -> Vector{T}`

Runs the anti-kt jet reconstruction algorithm but empties the given array of objects. See `anti_kt` for the non-mutating version.
"""
anti_kt!(objects::AbstractArray{T} where T; R=1, recombine=+) = sequential_jet_reconstruct!(objects, R=R, recombine=recombine)

"""
`anti_kt(objects; R=1, recombine=(x, y)->(x + y)) -> Vector`

Runs the anti-kt jet reconstruction algorithm. See `anti_kt!` for the mutating and thus slightly more efficient version.
"""
anti_kt(objects; R=1, recombine=+) = sequential_jet_reconstruct(objects, R=R, recombine=recombine)

end
