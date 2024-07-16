using Documenter
using DocumenterVitepress
using CairoMakie
using JetReconstruction

push!(LOAD_PATH, "../ext/")

include(joinpath(@__DIR__, "..", "ext", "JetVisualisation.jl"))

makedocs(sitename = "JetReconstruction.jl",
         format = MarkdownVitepress(repo = "github.com/JuliaHEP/JetReconstruction.jl"),
         pages = [
             "Home" => "index.md",
             "Examples" => "examples.md",
             "Reference" => Any["Public API" => "lib/public.md",
                                "Internal API" => "lib/internal.md",
                                "Visualisation API" => "lib/visualisation.md"]
         ])

deploydocs(repo = "github.com/JuliaHEP/JetReconstruction.jl.git",
           push_preview = true)
