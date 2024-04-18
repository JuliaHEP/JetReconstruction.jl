# Utility functions, which can be used by different top level scripts

"""Read HepMC3 event and keep final state particles (return PseudoJets)"""
function read_final_state_particles(fname; maxevents = -1, skipevents = 0)
	f = open(fname)

	events = Vector{PseudoJet}[]

	ipart = 1
	HepMC3.read_events(f, maxevents = maxevents, skipevents = skipevents) do parts
		input_particles = PseudoJet[]
		for p in parts
			if p.status == 1
				push!(
					input_particles,
					PseudoJet(p.momentum.x, p.momentum.y, p.momentum.z, p.momentum.t),
				)
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

"""Read HepMC3 event and keep final state particles (return LorentzVectors)"""
function read_final_state_particles_lv(fname; maxevents = -1, skipevents = 0)
	f = open(fname)

	events = Vector{LorentzVector{Float64}}[]

	ipart = 1
	HepMC3.read_events(f, maxevents = maxevents, skipevents = skipevents) do parts
		input_particles = LorentzVector{Float64}[]
		for p in parts
			if p.status == 1
				push!(
					input_particles,
					LorentzVector(p.momentum.t, p.momentum.x, p.momentum.y, p.momentum.z),
				)
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

function pseudojets2vectors(events::Vector{Vector{PseudoJet}})
	"""Convert events in PseudoJets into deep vectors"""
	event_vector = Vector{Vector{Vector{Float64}}}(undef, size(events)[1])
	for (ievent, event) in enumerate(events)
		jet_struct = Vector{Vector{Float64}}(undef, size(event)[1])
		for (ipart, initial_particle) in enumerate(event)
			jet_struct[ipart] = [
				initial_particle.px,
				initial_particle.py,
				initial_particle.pz,
				initial_particle.E,
			]
		end
		event_vector[ievent] = jet_struct
	end
	event_vector
end

"""
Return the list of jets passing a pt cut

The ptmin cut in these functions is slightly legacy as often the
input jets were already filtered on pt 
"""
function final_jets(jets::Vector{Vector{Float64}}, ptmin::AbstractFloat=0.0)
	count = 0
	final_jets = Vector{FinalJet}()
	sizehint!(final_jets, 6)
	for jet in jets
		dcut = ptmin^2
		p = PseudoJet(jet[1], jet[2], jet[3], jet[4])
		if p._pt2 > dcut
			count += 1
			push!(final_jets, FinalJet(rap(p), phi(p), sqrt(pt2(p))))
		end
	end
	final_jets
end

function final_jets(jets::Vector{PseudoJet}, ptmin::AbstractFloat=0.0)
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

function final_jets(jets::Vector{LorentzVector}, ptmin::AbstractFloat=0.0)
	count = 0
	final_jets = Vector{FinalJet}()
	sizehint!(final_jets, 6)
	dcut = ptmin^2
	for jet in jets
		if LorentzVectorHEP.pt(jet)^2 > dcut
			count += 1
			push!(final_jets, FinalJet(LorentzVectorHEP.rapidity(jet), LorentzVectorHEP.phi(jet), LorentzVectorHEP.pt(jet)))
		end
	end
	final_jets
end

function final_jets(jets::Vector{LorentzVectorCyl}, ptmin::AbstractFloat=0.0)
	count = 0
	final_jets = Vector{FinalJet}()
	sizehint!(final_jets, 6)
	dcut = ptmin^2
	for jet in jets
		if LorentzVectorHEP.pt(jet)^2 > dcut
			count += 1
			push!(final_jets, FinalJet(LorentzVectorHEP.eta(jet), LorentzVectorHEP.phi(jet), LorentzVectorHEP.pt(jet)))
		end
	end
	final_jets
end

## Faster findmin function
"""Find the lowest value in the array, returning the value and the index"""
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