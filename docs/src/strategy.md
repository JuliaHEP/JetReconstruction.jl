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

## Reusing N2Tiled storage

The ordinary `jet_reconstruct` and `tiled_jet_reconstruct` interfaces return an
independently owned `ClusterSequence`. This is the simplest interface to use,
but it does carry a performance penalty, which the workspace described below
overcomes.

High-throughput applications that process many events can avoid most per-event
temporary allocations with an `N2TiledWorkspace`:

```julia
workspace = N2TiledWorkspace()

for event in events
    with_n2tiled_reconstruction(
        workspace,
        event;
        algorithm = JetAlgorithm.AntiKt,
        R = 0.4,
    ) do clusterseq
        # This selection is independently owned and can outlive the callback.
        jets = inclusive_jets(clusterseq; ptmin = 5.0)

        # To retain the complete clustering sequence instead:
        # retained_clusterseq = deepcopy(clusterseq)
    end
end
```

The callback receives a complete `ClusterSequence`, so history, constituent,
and exclusive-jet queries remain available. Its jets and history borrow storage
from the workspace and are overwritten the next time that workspace is used.
You **must copy** values that will outlive the callback or use the ordinary
owning interface.

A workspace must not be shared concurrently or used reentrantly. Parallel
applications should create one workspace for each concurrent worker and keep
each workspace owned by that worker. `release_n2tiled_workspace_capacity!` can
be used after an unusually large event or when a long-lived worker should
release retained storage.

The [multithreaded N2Tiled example](https://github.com/JuliaHEP/JetReconstruction.jl/blob/main/examples/n2tiled-multithreaded.jl)
demonstrates dynamic event scheduling with one workspace owned by each
long-lived worker task.
