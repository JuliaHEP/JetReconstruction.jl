### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ cf88f2f9-bf19-47e2-ae04-b4476eb26efb
import Pkg

# ╔═╡ 24dd666c-49be-4c44-9b9d-d27128c1ffbb
Pkg.activate(".")

# ╔═╡ 0a1a8211-507b-4dc0-9060-0d18a5a695a4
using GLMakie

# ╔═╡ f7ec7400-bbb2-4f28-9f13-881dea775f41
using PlutoUI

# ╔═╡ 32b48d85-40cb-42c8-b3ad-ab613081da38
using JetReconstruction

# ╔═╡ dff6a188-2cbe-11ef-32d0-73c4c05efad2
md"""# Jet Reconstruction Visualisation

This Pluto script visualises the result of a jet reconstruction process. 

Use the sliders below to change the reconstructed event, the algorithm and the jet radius parameter."""

# ╔═╡ 0642cc6e-290e-4e01-9882-5c12d7585ead


# ╔═╡ 4e569f74-570b-4b30-9ea7-9cbc420f50f8
md"Event number:"

# ╔═╡ 83030910-8d00-4949-b69e-fa492b61db6b
@bind event_no PlutoUI.Slider(1:100, show_value=true)

# ╔═╡ 306f6e66-788c-4a95-a9df-3ac4c37cf776
md"``k_T`` algorithm power (-1=Anti-``k_T``, 0 = Cambridge/Aachen, 1=Inclusive kT)"

# ╔═╡ 864d7a31-b634-4edc-b6c4-0835fee53304
@bind power PlutoUI.Slider(-1:1, show_value=true)

# ╔═╡ 9f7988a5-041e-4a83-9da9-d48db34cea98
md"Jet radius parameter"

# ╔═╡ 187315d4-ac3c-4847-ae00-e7da1b27d63f
@bind radius PlutoUI.Slider(range(start=0.4, stop=2.0, step=0.1), show_value=true)

# ╔═╡ 79f24ec1-a63e-4e96-bd67-49661125be66
input_file = joinpath(dirname(pathof(JetReconstruction)), "..", "test", "data", "events.hepmc3.gz")

# ╔═╡ 7d7a8b11-19b3-4b83-a0b1-8201b74b588e
events::Vector{Vector{PseudoJet}} = read_final_state_particles(input_file, 
    maxevents = event_no, skipevents = event_no);

# ╔═╡ 2a899d67-71f3-4fe0-8104-7633a44a06a8
cs = jet_reconstruct(events[1], p=power, R=radius)

# ╔═╡ b5fd4e96-d073-4e5f-8de1-41addaa0dc3d
jetreco_vis = jetsplot(events[1], cs; Module = GLMakie)

# ╔═╡ 81b14eba-6981-4943-92fe-529c53b0ac35
display(jetreco_vis);

# ╔═╡ Cell order:
# ╟─dff6a188-2cbe-11ef-32d0-73c4c05efad2
# ╠═0642cc6e-290e-4e01-9882-5c12d7585ead
# ╠═cf88f2f9-bf19-47e2-ae04-b4476eb26efb
# ╠═24dd666c-49be-4c44-9b9d-d27128c1ffbb
# ╠═0a1a8211-507b-4dc0-9060-0d18a5a695a4
# ╠═f7ec7400-bbb2-4f28-9f13-881dea775f41
# ╠═32b48d85-40cb-42c8-b3ad-ab613081da38
# ╟─4e569f74-570b-4b30-9ea7-9cbc420f50f8
# ╟─83030910-8d00-4949-b69e-fa492b61db6b
# ╟─306f6e66-788c-4a95-a9df-3ac4c37cf776
# ╟─864d7a31-b634-4edc-b6c4-0835fee53304
# ╟─9f7988a5-041e-4a83-9da9-d48db34cea98
# ╟─187315d4-ac3c-4847-ae00-e7da1b27d63f
# ╟─79f24ec1-a63e-4e96-bd67-49661125be66
# ╠═7d7a8b11-19b3-4b83-a0b1-8201b74b588e
# ╠═2a899d67-71f3-4fe0-8104-7633a44a06a8
# ╠═b5fd4e96-d073-4e5f-8de1-41addaa0dc3d
# ╠═81b14eba-6981-4943-92fe-529c53b0ac35
