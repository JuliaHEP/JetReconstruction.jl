# Multithreaded Reconstruction

The ordinary `jet_reconstruct` and `tiled_jet_reconstruct` interfaces return an
independently owned `ClusterSequence`. They are the recommended interfaces for
ordinary, single-threaded reconstruction: their allocation overhead is small,
and the returned result can be retained freely.

## High-throughput multi-threaded reconstruction

High-throughput multi-threaded applications that reconstruct many independent
events in parallel can reuse N2Tiled storage with an `N2TiledWorkspace`. This
avoids most per-event temporary allocations, which greatly improves performance.
The cost of those allocations becomes more severe in multithreaded
applications.

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
from the workspace and are overwritten by the next reconstruction using that
workspace. Copy any data that must outlive the callback or use an ordinary
owning interface instead.

## Concurrent use

A workspace must not be shared concurrently or used reentrantly. Create one
workspace for each long-lived worker task and keep that workspace owned by the
task. The [multithreaded N2Tiled example](https://github.com/JuliaHEP/JetReconstruction.jl/blob/main/examples/n2tiled-multithreaded.jl)
shows dynamic scheduling of events in parallel while preserving that ownership.
