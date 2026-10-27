# Explicit ghost area reference

`jet-areas-fastjet.json.zst` was generated with FastJet 3.4.3. It stores the
physical inputs, exact generated ghost four-momenta, actual cell area, and final
jet momenta, scalar areas, and area four-vectors for AntiKt, CA, Kt, and GenKt
(`p = 1.5`), with E, pt, and pt2 recombination. The saved ghost inputs make the
reference independent of Julia's RNG implementation. The event includes a
physical pair across the azimuthal seam and a jet near the ghost coverage edge.

Regenerate from the repository root with FastJet installed:

```sh
c++ -std=c++11 $(fastjet-config --cxxflags) test/data/jet-areas-fastjet.cc \
    $(fastjet-config --libs) -o /tmp/jet-areas-fastjet
/tmp/jet-areas-fastjet | zstd -o test/data/jet-areas-fastjet.json.zst -f
```

The generator uses FastJet's ghost generator with random status `[12345, 67890]`
and passes the same explicit ghosts to every reconstruction. Tests use the saved
inputs with both Julia strategies and compare against the FastJet outputs.
