# Utility functions, which can be used by different top level scripts

using CodecZlib
using CodecZstd

"""
    open_with_stream(fname::AbstractString)

Open a file with a stream decompressor if it is compressed with gzip or zstd,
otherwise as a normal file.
"""
open_with_stream(fname::AbstractString) = begin
    if endswith(fname, ".gz")
        f = GzipDecompressorStream(open(fname))
    elseif endswith(fname, ".zst")
        f = ZstdDecompressorStream(open(fname))
    else
        f = open(fname)
    end
    f
end

"""
    read_final_state_particles(fname; maxevents = -1, skipevents = 0, T=PseudoJet)

Reads final state particles from a file and returns them as a vector of type T.

# Arguments
- `fname`: The name of the HepMC3 ASCII file to read particles from. If the file
  is gzipped, the function will automatically decompress it.
- `maxevents=-1`: The maximum number of events to read. -1 means all events will
  be read.
- `skipevents=0`: The number of events to skip before an event is included.
- `T=PseudoJet`: The type of object to construct and return.

# Returns
A vector of vectors of T objects, where each inner vector represents all
the particles of a particular event. In particular T can be `PseudoJet` or
a `LorentzVector` type. Note, if T is not `PseudoJet`, the order of the
arguments in the constructor must be `(t, x, y, z)`.
"""
function read_final_state_particles(fname; maxevents = -1, skipevents = 0, T = PseudoJet)
    f = open_with_stream(fname)
    events = Vector{T}[]

    ipart = 1
    HepMC3.read_events(f, maxevents = maxevents, skipevents = skipevents) do parts
        input_particles = T[]
        for p in parts
            if p.status == 1
                # Annoyingly PseudoJet and LorentzVector constructors
                # disagree on the order of arguments...
                if T <: FourMomentum
                    particle = T(p.momentum.x, p.momentum.y, p.momentum.z, p.momentum.t)
                else
                    particle = T(p.momentum.t, p.momentum.x, p.momentum.y, p.momentum.z)
                end
                push!(input_particles, particle)
            end
        end
        push!(events, input_particles)
        ipart += 1
    end
    close(f)

    @info "Total Events: $(length(events))"
    @debug events
    events
end

"""
    final_jets(jets::Vector{PseudoJet}, ptmin::AbstractFloat=0.0)

This function takes a vector of `PseudoJet` objects and a minimum transverse
momentum `ptmin` as input. It returns a vector of `FinalJet` objects that
satisfy the transverse momentum condition.

# Arguments
- `jets::Vector{PseudoJet}`: A vector of `PseudoJet` objects representing the
  input jets.
- `ptmin::AbstractFloat=0.0`: The minimum transverse momentum required for a jet
  to be included in the final jets vector.

# Returns
A vector of `FinalJet` objects that satisfy the transverse momentum condition.
"""
function final_jets(jets::Vector{PseudoJet}, ptmin::AbstractFloat = 0.0)
    count = 0
    final_jets = Vector{FinalJet}()
    sizehint!(final_jets, 6)
    for jet in jets
        dcut = ptmin^2
        if pt2(jet) > dcut
            count += 1
            push!(final_jets, FinalJet(rapidity(jet), phi(jet), sqrt(pt2(jet))))
        end
    end
    final_jets
end

"""Specialisation for final jets from LorentzVectors (TODO: merge into more general function)"""
function final_jets(jets::Vector{LorentzVector}, ptmin::AbstractFloat = 0.0)
    count = 0
    final_jets = Vector{FinalJet}()
    sizehint!(final_jets, 6)
    dcut = ptmin^2
    for jet in jets
        if LorentzVectorHEP.pt(jet)^2 > dcut
            count += 1
            push!(final_jets,
                  FinalJet(LorentzVectorHEP.rapidity(jet), LorentzVectorHEP.phi(jet),
                           LorentzVectorHEP.pt(jet)))
        end
    end
    final_jets
end

"""Specialisation for final jets from LorentzVectorCyl (TODO: merge into more general function)"""
function final_jets(jets::Vector{LorentzVectorCyl}, ptmin::AbstractFloat = 0.0)
    count = 0
    final_jets = Vector{FinalJet}()
    sizehint!(final_jets, 6)
    dcut = ptmin^2
    for jet in jets
        if LorentzVectorHEP.pt(jet)^2 > dcut
            count += 1
            push!(final_jets,
                  FinalJet(LorentzVectorHEP.eta(jet), LorentzVectorHEP.phi(jet),
                           LorentzVectorHEP.pt(jet)))
        end
    end
    final_jets
end

"""
    fast_findmin(dij, n)

Find the minimum value and its index in the first `n` elements of the `dij`
array. The use of `@turbo` macro gives a significant performance boost.

# Arguments
- `dij`: An array of values.
- `n`: The number of elements to consider in the `dij` array.

# Returns
- `dij_min`: The minimum value in the first `n` elements of the `dij` array.
- `best`: The index of the minimum value in the `dij` array.
"""
function fast_findmin end

if Sys.ARCH == :aarch64
    fast_findmin(dij, n) = _naive_fast_findmin(@view(dij[begin:n]))
else
    function fast_findmin(dij, n)
        if n <= 8
            return _naive_fast_findmin(@view(dij[begin:n]))
        else
            return _simd_fast_findmin(dij, n)
        end
    end
end

function _naive_fast_findmin(dij)
    x = @fastmath foldl(min, dij)
    i = findfirst(==(x), dij)::Int
    x, i
end

function _simd_fast_findmin(dij::DenseVector{T}, n) where {T}
    laneIndices = SIMD.Vec{8, Int}((1, 2, 3, 4, 5, 6, 7, 8))
    minvals = SIMD.Vec{8, T}(Inf)
    min_indices = SIMD.Vec{8, Int}(0)

    n_batches, remainder = divrem(n, 8)
    lane = VecRange{8}(0)
    i = 1
    @inbounds @fastmath for _ in 1:n_batches
        dijs = dij[lane + i]
        predicate = dijs < minvals
        minvals = vifelse(predicate, dijs, minvals)
        min_indices = vifelse(predicate, laneIndices, min_indices)

        i += 8
        laneIndices += 8
    end

    # last batch
    back_track = 8 - remainder
    i -= back_track
    laneIndices -= back_track

    dijs = dij[lane + i]
    predicate = dijs < minvals
    minvals = vifelse(predicate, dijs, minvals)
    min_indices = vifelse(predicate, laneIndices, min_indices)

    min_value = SIMD.minimum(minvals)
    min_index = @inbounds min_value == minvals[1] ? min_indices[1] :
                          min_value == minvals[2] ? min_indices[2] :
                          min_value == minvals[3] ? min_indices[3] :
                          min_value == minvals[4] ? min_indices[4] :
                          min_value == minvals[5] ? min_indices[5] :
                          min_value == minvals[6] ? min_indices[6] :
                          min_value == minvals[7] ? min_indices[7] : min_indices[8]

    return min_value, min_index
end
