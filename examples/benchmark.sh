#! /bin/sh
#
# Quick and dirty set of benchmarks for the most important cases
echo "pp 14TeV Tiled"
julia --project instrumented-jetreco.jl --algorithm=AntiKt -R 0.4 ../test/data/events.pp13TeV.hepmc3.gz -S N2Tiled -m 16

echo "pp 14 TeV Plain"
julia --project instrumented-jetreco.jl --algorithm=AntiKt -R 0.4 ../test/data/events.pp13TeV.hepmc3.gz -S N2Plain -m 16

echo "ee H Durham"
julia --project instrumented-jetreco.jl --algorithm=Durham ../test/data/events.eeH.hepmc3.gz -m 16
