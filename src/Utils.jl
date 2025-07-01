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
    read_final_state_particles(fname, ::Type{T} = PseudoJet; maxevents = -1, skipevents = 0) where {T}

Reads final state particles from a file and returns them as a vector of type T.

# Arguments
- `fname`: The name of the HepMC3 ASCII file to read particles from. If the file
  is gzipped, the function will automatically decompress it.
- `::Type{T}=PseudoJet`: The type of object to construct and return.
- `maxevents=-1`: The maximum number of events to read. -1 means all events will
  be read.
- `skipevents=0`: The number of events to skip before an event is included.

# Returns
A vector of vectors of T objects, where each inner vector represents all
the particles of a particular event. In particular T can be `PseudoJet` or
a `LorentzVector` type. Note, if T is not `PseudoJet`, the order of the
arguments in the constructor must be `(t, x, y, z)`.
"""
function read_final_state_particles(fname, ::Type{T} = PseudoJet; maxevents = -1,
                                    skipevents = 0) where {T}
    f = open_with_stream(fname)
    events = Vector{T}[]

    ipart = 1
    HepMC3.read_events(f, maxevents = maxevents, skipevents = skipevents) do parts
        input_particles = T[]
        particle_index = 1
        for p in parts
            if p.status == 1
                # Annoyingly PseudoJet and LorentzVector constructors
                # disagree on the order of arguments...
                if T <: FourMomentum
                    particle = T(p.momentum.x, p.momentum.y, p.momentum.z, p.momentum.t;
                                 cluster_hist_index = particle_index)
                else
                    particle = T(p.momentum.t, p.momentum.x, p.momentum.y, p.momentum.z)
                end
                push!(input_particles, particle)
                particle_index += 1
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
    final_jets(jets::Vector{T}, ptmin::AbstractFloat=0.0)

This function takes a vector of `T` objects, which should be jets, and a minimum
transverse momentum `ptmin` as input. It returns a vector of `FinalJet` objects
that satisfy the transverse momentum condition.

# Arguments
- `jets::Vector{T}`: A vector of `T` objects representing the input jets.
- `ptmin::AbstractFloat=0.0`: The minimum transverse momentum required for a jet
  to be included in the final jets vector.

# Returns
A vector of `FinalJet` objects that satisfy the transverse momentum condition.
"""
function final_jets(jets::Vector{T}, ptmin::AbstractFloat = 0.0) where {T}
    count = 0
    final_jets = Vector{FinalJet}()
    dcut = ptmin^2
    for jet in jets
        if pt2(jet) > dcut
            count += 1
            push!(final_jets, FinalJet(jet))
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
fast_findmin(dij, n) = begin
    # findmin(@inbounds @view dij[1:n])
    best = 1
    @inbounds dij_min = dij[1]
    @turbo for here in 2:n
        newmin = dij[here] < dij_min
        best = newmin ? here : best
        dij_min = newmin ? dij[here] : dij_min
    end
    dij_min, best
end
