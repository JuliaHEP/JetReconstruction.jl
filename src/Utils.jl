# Utility functions, which can be used by different top level scripts

using CodecZlib

"""
    read_final_state_particles(fname; maxevents = -1, skipevents = 0, T=PseudoJet)

Reads final state particles from a file and returns them as a vector of type T.

# Arguments
- `fname`: The name of the HepMC3 ASCII file to read particles from. If the file
  is gzipped, the function will automatically decompress it.
- `maxevents=-1`: The maximum number of events to read. -1 means all events will
  be read.
- `skipevents=0`: The number of events to skip before an event is included.
- `T=PseudoJet`: The type of object to contruct and return.

# Returns
A vector of vectors of T objects, where each inner vector represents all
the particles of a particular event. In particular T can be `PseudoJet` or
a `LorentzVector` type. Note, if T is not `PseudoJet`, the order of the
arguments in the constructor must be `(t, x, y, z)`.
"""
function read_final_state_particles(fname; maxevents = -1, skipevents = 0, T = PseudoJet)
    f = open(fname)
    if endswith(fname, ".gz")
        @debug "Reading gzipped file $fname"
        f = GzipDecompressorStream(f)
    end

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
array. The use of `@turbo` macro gives a significiant performance boost.

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
