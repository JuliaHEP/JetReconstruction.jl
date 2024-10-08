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
             "EDM4hep" => "EDM4hep.md",
             "Visualisation" => "visualisation.md",
             "Reference" => Any["Public API" => "lib/public.md",
                                "Internal API" => "lib/internal.md"]
         ])

deploydocs(repo = "github.com/JuliaHEP/JetReconstruction.jl.git",
           devbranch = "main",
           devurl = "dev",
           versions = ["stable" => "v^", "v#.#", "dev" => "dev"],
           push_preview = true)
