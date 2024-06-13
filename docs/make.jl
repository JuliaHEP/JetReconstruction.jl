using Documenter
using JetReconstruction

makedocs(sitename = "JetReconstruction.jl")

deploydocs(repo = "github.com/JuliaHEP/JetReconstruction.jl.git",
           push_preview = true)
