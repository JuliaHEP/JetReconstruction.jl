# Algorithm Strategy

For the ``pp`` algorithms three strategies are available for the different
algorithms, which can be specified by passing the named argument `strategy=...`
to the reconstruction.

| Strategy Name | Notes | Interface |
|---|---|---|
| `RecoStrategy.Best` | Dynamically switch strategy based on input particle density | `jet_reconstruct` |
| `RecoStrategy.N2Plain` | Global matching of particles at each interaction (works well for low $N$) | `plain_jet_reconstruct` |
| `RecoStrategy.N2Tiled` | Use tiles of radius $R$ to limit search space (works well for higher $N$) | `tiled_jet_reconstruct` |

Generally one can use the `jet_reconstruct` interface, shown above, as the
*Best* strategy safely as the overhead is extremely low. That interface supports
a `strategy` option to switch to a different option.

For ``e^+e^-`` algorithms particle densities are low, so the only
implementation for these algorithms is effectively of the same type as
`N2Plain`.
