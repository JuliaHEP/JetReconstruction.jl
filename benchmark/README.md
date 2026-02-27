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
julia --project=benchmark benchmark/runbenchmarks.jl --verbose benchmark_results.json
```

To control how many threads the benchmark process uses, pass `--threads N`.
This is forwarded via `BenchmarkConfig` to the subprocess that actually runs
the benchmarks:

```sh
julia --project=benchmark benchmark/runbenchmarks.jl --threads 4 --verbose benchmark_results.json
```
