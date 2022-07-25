# JetReconstruction

[![Build Status](https://github.com/gojakuch/JetReconstruction.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/gojakuch/JetReconstruction.jl/actions/workflows/CI.yml?query=branch%3Amain)

This code is not a registered Julia package yet and it is not optimised. For now, using it is not recommended. If you want to use it, the valid option would be to change the `main.jl` file how you want and run it. That file should either start with
```julia
using Revise; import Pkg # you need to have Revise installed
Pkg.activate(".")
using JetReconstruction
```
or with
```julia
include("src/JetReconstruction.jl")
using .JetReconstruction
```
There is some documentation provided for functions and submodules. Once everything is working properly and efficiently, this `README.md` will contain more details on usage and the module might get registered.

## Plotting
![illustration](img/illustration.jpeg)

To visualise the clustered jets as a 3d bar plot (see illustration above) the `PyPlot` package is used (for now), as native Julia plotting libraries do not provide easy access to such functionality. Use the `jetsplot` function. See its documentation for more.
