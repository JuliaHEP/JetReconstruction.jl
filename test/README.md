# Tests for JetReconstruction

This directory contains the unit tests for the JetReconstruction package.

To run all of the tests do:

```sh
julia --project test/runtests.jl
```

or

```sh
julia --project
julia> ]
(JetReconstruction) pkg> test
...
```

Tests are factorised into individual files as `test-*.jl`. By use of `common.jl` all of these tests can be run standalone, provided that `TestEnv.jl` is installed in the default environment, e.g.,

```sh
julia --project test/test-pp-reconstruction.jl
```
