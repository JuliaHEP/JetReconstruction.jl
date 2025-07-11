using Documenter
using CairoMakie
using EDM4hep
using JetReconstruction

push!(LOAD_PATH, "../ext/")

include(joinpath(@__DIR__, "..", "ext", "JetVisualisation.jl"))
include(joinpath(@__DIR__, "..", "ext", "EDM4hepJets.jl"))

makedocs(sitename = "JetReconstruction.jl",
         clean = false,
         pages = [
             "Home" => "index.md",
             "Examples" => "examples.md",
             "Particle Inputs" => "particles.md",
             "Reconstruction Strategies" => "strategy.md",
             "Substructure" => "substructure.md",
             "Jet Helpers" => "helpers.md",
             "EDM4hep" => "EDM4hep.md",
             "Recombination Schemes" => "recombination.md",
             "SoftKiller" => "softkiller.md",
             "Visualisation" => "visualisation.md",
             "Contributing" => "contributing.md",
             "Reference Docs" => Any["Public API" => "lib/public.md",
                                     "Internal API" => "lib/internal.md"]
         ])

deploydocs(repo = "github.com/JuliaHEP/JetReconstruction.jl.git",
           devbranch = "main",
           devurl = "dev",
           versions = ["stable" => "v^", "v#.#", "dev" => "dev"],
           push_preview = true)
