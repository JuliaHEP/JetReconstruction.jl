#! /bin/sh
# Simple test file to check that examples are running correctly
echo "pp visualisation"
julia --project visualise-jets.jl -A AntiKt -R 1.0 ../../test/data/events.pp13TeV.hepmc3.gz test-pp.png

echo "e+e- visualisation"
julia --project visualise-jets.jl -A Durham ../../test/data/events.eeH.hepmc3.gz test-ee.png

echo "pp animation"
julia --project animate-reconstruction.jl -A AntiKt -R 1.0 ../../test/data/events.pp13TeV.hepmc3.gz test-pp.mp4

echo "e+e- animation"
julia --project animate-reconstruction.jl -A Durham ../../test/data/events.eeH.hepmc3.gz test-eeH.mp4
