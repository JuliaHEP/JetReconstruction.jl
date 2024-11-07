# Generating and Previewing Documentation

## Adding Documentation

`Documenter.jl` is used to generate package documentation. 

The general guidelines are to have documentation split into meaningful sections,
none of which are too long. If a new functionality is added then create a new
documentation page in `src`. Then map that markdown file to the a section in the
documentation in `make.jl` in the `pages` map.

## Previewing Documenation

From the `docs` directory:

```julia
julia --project
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.11.1 (2024-10-16)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> include("make.jl")
[ Info: SetupBuildDirectory: setting up build directory.
[ Info: Doctest: running doctests.
[ Info: ExpandTemplates: expanding markdown templates.
[ Info: CrossReferences: building cross-references.
[ Info: CheckDocument: running document checks.
[ Info: Populate: populating indices.
[ Info: RenderDocument: rendering document.
[ Info: HTMLWriter: rendering HTML pages.
┌ Warning: Generated HTML over size_threshold_warn limit: lib/internal.md
│     Generated file size: 128.16 (KiB)
│     size_threshold_warn: 100.0 (KiB)
│     size_threshold:      200.0 (KiB)
│     HTML file:           lib/internal/index.html
└ @ Documenter.HTMLWriter ~/.julia/packages/Documenter/C1XEF/src/html/HTMLWriter.jl:1828
[ Info: Automatic `version="0.4.3"` for inventory from ../Project.toml
┌ Warning: Documenter could not auto-detect the building environment. Skipping deployment.
└ @ Documenter ~/.julia/packages/Documenter/C1XEF/src/deployconfig.jl:76

julia> using LiveServer

julia> LiveServer.serve(dir="build/")
✓ LiveServer listening on http://localhost:8000/ ...
  (use CTRL+C to shut down)
```
