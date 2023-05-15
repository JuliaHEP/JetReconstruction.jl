# Utility functions, which can be used by different top level scripts

"""Read HepMC3 event and keep final state particles"""
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

	@info "Total Events: $(length(events))"
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

function final_jets(jets::Vector{Vector{Float64}}, ptmin::AbstractFloat)
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
