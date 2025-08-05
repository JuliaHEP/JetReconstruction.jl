#! /bin/sh

# Plain softkiller example
julia --project softkiller_runtime.jl --pileup-maxevents=100 --eventno=4 --grid-size=0.4 --algorithm=Kt ../../test/data/sk_example_HS.hepmc.zst ../../test/data/sk_example_PU.hepmc.zst

# Softkiller with graphics
julia --project softkiller_plots.jl --pileup-maxevents=100 --eventno=4 --grid-size=0.4 --algorithm=Kt ../../test/data/sk_example_HS.hepmc.zst ../../test/data/sk_example_PU.hepmc.zst
