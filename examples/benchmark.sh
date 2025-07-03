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
