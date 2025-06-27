#! /bin/sh
# Simple test file to check that examples are running correctly
echo "pp visualisation"
julia --project visualise-jets.jl -A AntiKt -R 1.0 ../../test/data/events.pp13TeV.hepmc3.zst test-pp.png

echo "e+e- visualisation"
julia --project visualise-jets.jl -A Durham ../../test/data/events.eeH.hepmc3.zst test-ee.png

echo "pp animation"
julia --project animate-reconstruction.jl -A AntiKt -R 1.0 ../../test/data/events.pp13TeV.hepmc3.zst test-pp.mp4

echo "e+e- animation"
julia --project animate-reconstruction.jl -A Durham ../../test/data/events.eeH.hepmc3.zst test-eeH.mp4
