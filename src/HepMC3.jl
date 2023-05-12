# Copyright (c) 2022, Philippe Gras
#
#----------------------------------------------------------------------
# This file is part of AntiKt.jl.
#
#  AntiKt.jl is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  The algorithms that underlie FastJet have required considerable
#  development. They are described in the original FastJet paper,
#  hep-ph/0512210 and in the manual, arXiv:1111.6097. If you use
#  FastJet as part of work towards a scientific publication, please
#  quote the version you use and include a citation to the manual and
#  optionally also to hep-ph/0512210.
#
#  AntiKet.jl is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with FastJet. If not, see <http:#www.gnu.org/licenses/>.
#----------------------------------------------------------------------

""" Module providing a minimal support to read HepMC3 ascii files.
"""
module HepMC3

struct Momentum4{T}
    px::T
    py::T
    pz::T
    e::T
end

struct Particle{T}
    momentum::Momentum4{T}
    status::Integer
    pdgid::Integer
    barcode::Integer
    vertex::Integer
end

Particle{T}() where T = Particle(Momentum4{T}(0., 0., 0., 0.), 0, 0, 0, 0)

""" Read a [HepMC3](https://doi.org/10.1016/j.cpc.2020.107310) ascii file.

    Each event is passed to the provided function f as a vector of Particles. A
    maximum number of events to read (value -1 to read all availble events) and
    a number of events to skip at the beginning of the file can be provided.
"""
function read_events(f, fin; maxevents=-1, skipevents=0)
    T = Float64
    particles = Particle{T}[]
    ievent = 0
    ipart = 0
    toskip = skipevents
    
    for (il, l) in enumerate(eachline(fin))
        if occursin(r"HepMC::.*-END_EVENT_LISTING", l)
            break
        end

        tok = split(l)

        if tok[1] == "E"
            ievent += 1
            (maxevents >= 0 && ievent > maxevents) && break
            if ievent > 1 && toskip == 0
                f(particles)
            end
            if toskip > 0
                toskip -= 1
            end
            resize!(particles, 0)
            ipart = 0
        elseif tok[1] == "P" && toskip == 0
            ipart += 1
            if ipart > length(particles)
                push!(particles, Particle{T}())
            end
            barcode = parse(Int, tok[2])
            vertex = parse(Int, tok[3])
            pdgid = parse(Int, tok[4])
            px = parse(T, tok[5])
            py = parse(T, tok[6])
            pz = parse(T, tok[7])
            e =  parse(T, tok[8])
            status = parse(Int, tok[10])
            push!(particles, Particle{T}(Momentum4(px,py,pz,e), status, pdgid, barcode, vertex))
        end
    end
    #processing the last event:
    ievent > 0 && f(particles)
end
end
