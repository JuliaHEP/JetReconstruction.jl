#! /bin/sh
#
# Quick and dirty set of benchmarks for the most important cases
trials=${1:-16}

echo "pp 14TeV Tiled"
julia --project instrumented-jetreco.jl --algorithm=AntiKt -R 0.4 ../test/data/events.pp13TeV.hepmc3.zst -S N2Tiled -m $trials

echo "pp 14 TeV Plain"
julia --project instrumented-jetreco.jl --algorithm=AntiKt -R 0.4 ../test/data/events.pp13TeV.hepmc3.zst -S N2Plain -m $trials

echo "ee H Durham"
julia --project instrumented-jetreco.jl --algorithm=Durham ../test/data/events.eeH.hepmc3.zst -m $trials

echo "ee Valencia"
julia --project ./instrumented-jetreco.jl --algorithm=Valencia --gamma=1.2 --power=1.2 -R 0.8 ../test/data/events.eeH.hepmc3.zst -m $trials
