# Mini-benchmarks for JetReconstruction

This directory contains mini-benchmarks for the `JetReconstruction.jl` package. A more comprehensive suite of benchmarks can be found in [JetReconstructionBenchmarks.jl](https://github.com/graeme-a-stewart/JetReconstructionBenchmarks.jl).

## PkgBenchmark-style benchmarks

The [benchmark/benchmarks.jl](benchmark/benchmarks.jl) benchmarks follow the standard [`PkgBenchmark.jl`](https://juliaci.github.io/PkgBenchmark.jl/stable/) structure and can be run locally with your favourite `PkgBenchmark.jl`-compatible tool, for example:

```julia 
using PkgBenchmark
import JetReconstruction
benchmarkpkg(JetReconstruction)
```

Alternatively, the benchmarks can be run directly from the command line with:

```sh
julia -t auto --project=benchmark benchmark/runbenchmarks.jl --verbose benchmark_results.json
```

The `-t auto` flag tells Julia to use all available CPU threads. The runner
reports the active thread count at startup so results are always self-describing.
