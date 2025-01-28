# Tests for JetReconstruction

This directory contains the unit tests for the JetReconstruction package.

As there are dependencies that are only used for tests, there is a specific
`Project.toml` file for testing. However, **make sure that the
`JetReconstruction` package is in develop mode**.

To run all of the tests do:

```sh
julia --project=test runtests.jl
```

or

```sh
julia --project
julia> ]
(JetReconstruction) pkg> test
...
```

Tests are factorised into individual files as `test-*.jl`. By use of `common.jl` all of these tests can be run standalone, e.g.,

```sh
julia --project=test test-pp-reconstruction.jl
```
