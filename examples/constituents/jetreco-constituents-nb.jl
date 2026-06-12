### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 21688366-1c52-4c12-a4b7-8db2375bd368
using Pkg

# ╔═╡ 727c300d-ad3d-437d-bdeb-88b18cdcee04
Pkg.activate(".")

# ╔═╡ 32b48d85-40cb-42c8-b3ad-ab613081da38
using JetReconstruction

# ╔═╡ dff6a188-2cbe-11ef-32d0-73c4c05efad2
md"""# Jet Reconstruction Constituents Example

Perform a simple reconstruction example and show how to retrieve constituent jets.

Note that the *local* package setup is used here, to get the current version of `JetReconstruction`."""

# ╔═╡ 79f24ec1-a63e-4e96-bd67-49661125be66
input_file = joinpath(dirname(pathof(JetReconstruction)), "..", "test", "data",
                      "events.pp13TeV.hepmc3.zst")

# ╔═╡ 7d7a8b11-19b3-4b83-a0b1-8201b74b588e
events = read_final_state_particles(input_file)

# ╔═╡ fe05d243-6d19-47cc-814b-0f44e5424233
md"Event to pick"

# ╔═╡ 6bbf7e6a-30a2-4917-bc36-e0da4f2cd098
event_no = 1

# ╔═╡ 2a899d67-71f3-4fe0-8104-7633a44a06a8
cluster_seq = jet_reconstruct(events[event_no]; algorithm = JetAlgorithm.Kt, R = 1.0)

# ╔═╡ 043cc484-e537-409c-8aa2-e4904b5dc283
md"Retrieve the exclusive pj_jets, but as `PseudoJet` types"

# ╔═╡ 0d8d4664-915f-4f28-9d5a-6e03cb8d7d8b
pj_jets = inclusive_jets(cluster_seq, PseudoJet; ptmin = 5.0)

# ╔═╡ 0bd764f9-d427-43fc-8342-603b6759ec8f
md"Get the constituents of the first jet"

# ╔═╡ 46a64c6f-51d7-4083-a953-ecc76882f21e
my_constituents = JetReconstruction.constituents(pj_jets[1], cluster_seq)

# ╔═╡ 300879ca-b53d-40b3-864a-1d46f2094123
begin
    println("Constituents of jet number $(event_no):")
    for c in my_constituents
        println(" $c")
    end
end

# ╔═╡ c9ce9c76-82ef-42ff-bb2e-3b3b8085d8bc
begin
	my_constituent_indexes = JetReconstruction.constituent_indexes(pj_jets[1], cluster_seq)
	println("\nConsitituent indexes for jet number $(event_no): $my_constituent_indexes")
	for i in my_constituent_indexes
	    println("  Constituent jet $i: $(events[1][i])")
	end
end

# ╔═╡ Cell order:
# ╟─dff6a188-2cbe-11ef-32d0-73c4c05efad2
# ╠═21688366-1c52-4c12-a4b7-8db2375bd368
# ╠═727c300d-ad3d-437d-bdeb-88b18cdcee04
# ╠═32b48d85-40cb-42c8-b3ad-ab613081da38
# ╠═79f24ec1-a63e-4e96-bd67-49661125be66
# ╠═7d7a8b11-19b3-4b83-a0b1-8201b74b588e
# ╟─fe05d243-6d19-47cc-814b-0f44e5424233
# ╠═6bbf7e6a-30a2-4917-bc36-e0da4f2cd098
# ╠═2a899d67-71f3-4fe0-8104-7633a44a06a8
# ╟─043cc484-e537-409c-8aa2-e4904b5dc283
# ╠═0d8d4664-915f-4f28-9d5a-6e03cb8d7d8b
# ╟─0bd764f9-d427-43fc-8342-603b6759ec8f
# ╠═46a64c6f-51d7-4083-a953-ecc76882f21e
# ╠═300879ca-b53d-40b3-864a-1d46f2094123
# ╠═c9ce9c76-82ef-42ff-bb2e-3b3b8085d8bc
