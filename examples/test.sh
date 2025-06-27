#! /bin/sh
# Simple test file to check that examples are running correctly
echo "pp 14TeV Tiled benchmark"
julia --project instrumented-jetreco.jl --algorithm=AntiKt -R 0.4 ../test/data/events.pp13TeV.hepmc3.zst -S N2Tiled -m 16

echo "pp 14 TeV Plain profile"
julia --project instrumented-jetreco.jl --algorithm=AntiKt -R 0.4 ../test/data/events.pp13TeV.hepmc3.zst -S N2Plain -m 2 --profile test

echo "pp 14TeV CA"
julia --project instrumented-jetreco.jl --algorithm=CA -R 1.0 ../test/data/events.pp13TeV.hepmc3.zst

echo "pp 14TeV Kt"
julia --project instrumented-jetreco.jl --algorithm=Kt -R 1.0 ../test/data/events.pp13TeV.hepmc3.zst

echo "ee H Durham allocation test"
julia --project instrumented-jetreco.jl --algorithm=Durham --alloc ../test/data/events.eeH.hepmc3.zst

echo "pp basic reco test"
julia --project jetreco.jl --algorithm=CA -R 1.0 ../test/data/events.pp13TeV.hepmc3.zst

echo "ee basic reco test"
julia --project jetreco.jl --algorithm=EEKt -R 1.0 -p -1.0 ../test/data/events.eeH.hepmc3.zst

