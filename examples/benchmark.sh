#! /bin/sh
#
# Quick and dirty set of benchmarks for the most important cases
trials=${1:-16}
numtype=${2:-Float64}

bold="\033[1m"
boldoff="\033[0m"

echo "${bold}pp 14TeV Tiled for ${numtype}${boldoff}"
julia --project instrumented-jetreco.jl --algorithm=AntiKt -R 0.4 ../test/data/events.pp13TeV.hepmc3.zst -S N2Tiled -m $trials --numtype $numtype
echo

echo "${bold}pp 14 TeV Plain for ${numtype}${boldoff}"
julia --project instrumented-jetreco.jl --algorithm=AntiKt -R 0.4 ../test/data/events.pp13TeV.hepmc3.zst -S N2Plain -m $trials --numtype $numtype
echo

echo "${bold}ee H Durham for ${numtype}${boldoff}"
julia --project instrumented-jetreco.jl --algorithm=Durham ../test/data/events.eeH.hepmc3.zst -m $trials --numtype $numtype
echo

echo "${bold}ee Valencia for ${numtype}${boldoff}"
julia --project ./instrumented-jetreco.jl --algorithm=Valencia --gamma=1.2 --power=1.2 -R 0.8 ../test/data/events.eeH.hepmc3.zst -m $trials --numtype $numtype
