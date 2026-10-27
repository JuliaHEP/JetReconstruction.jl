# Active Jet Areas

Jet areas measure the susceptibility of a jet to diffuse soft radiation.
`jet_reconstruct_area` adds very soft, massless **ghosts** over a rapidity–azimuth
grid and follows their contribution through clustering. The scalar area is the
number of ghost constituents times their actual cell area. Area four-vectors
are propagated using the selected recombination scheme.

```julia
using JetReconstruction
using Random

particles = [PseudoJet(; pt = 100.0, rap = 0.0, phi = 0.0)]
spec = GhostedAreaSpec(ghost_maxrap = 1.5, ghost_area = 0.01)
csa = jet_reconstruct_area(particles; algorithm = JetAlgorithm.AntiKt,
                           R = 0.4, area_spec = spec, rng = MersenneTwister(12345))
jets = inclusive_jets(csa; ptmin = 1.0)
jet = only(jets)
area(jet, csa)
area_4vector(jet, csa)
```

For an isolated anti-kt jet away from the acceptance boundary, the scalar area
is approximately `π * R^2`. All pp algorithms and reconstruction strategies are
supported. `GenKt` also requires `p`. `RecoStrategy.Best` is the default and
normally selects tiled clustering for the ghost multiplicity. Physical inputs
are copied, so their history indices need not be assigned beforehand.

## Coverage and Precision

Set `ghost_maxrap` at least `R` beyond the rapidities of the jets of interest.
Jets near the ghost coverage boundary have truncated areas. Reducing
`ghost_area` improves spatial resolution at the cost of more ghosts and work.
The actual cell area is adjusted so the grid fits its rapidity range and `2π`.
`total_area(csa)` is the full grid area, `4π * ghost_maxrap`.

The default ghosts have mean transverse momentum `1e-100`. Check
`has_dangerous_particles(csa)` for inputs whose momenta may be affected by the
ghosts; this uses the same `eps(Float64)^2` softness criterion as FastJet.
`max_ghost_pt2(csa)` reports the largest initial ghost momentum squared.

Use `rng` for control of the random stream. Recreating an RNG with the same seed
reproduces the distribution; reusing the same RNG advances its state. The
`GhostedAreaSpec(seed = ...)` convenience option creates a local RNG per call,
and an explicit `rng` overrides it.

## Jets, Ghosts, and History

`inclusive_jets(csa)` returns `PseudoJet`s and filters pure ghost jets by default.
Use `include_pure_ghosts = true` to include them, and `is_pure_ghost(jet, csa)`
to identify them. The sum of scalar areas over all inclusive jets, including
pure ghost jets and without a momentum cut, equals `total_area(csa)`.

Area queries require a `PseudoJet` from the same area sequence. Conversions such
as `inclusive_jets(csa, LorentzVector)` discard its history index. The underlying
sequence is accessible with `cluster_sequence(csa)` for history and constituent
queries. Constituent indices beyond `n_hard_particles(csa)` identify ghosts;
the remaining indices address the original physical inputs.

## Scope

This implements active areas with **explicit ghosts**, using a single ghost
distribution (`repeat = 1`). `area_error` returns zero because no uncertainty
from repeated distributions is estimated. Passive areas, Voronoi areas,
background density estimation, and pile-up subtraction are outside this API.

The grid and area propagation follow [FastJet's explicit-ghost area method](https://fastjet.fr/repo/fastjet-doc-3.4.3.pdf).
