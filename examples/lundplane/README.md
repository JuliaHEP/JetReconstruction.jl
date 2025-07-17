# Examples for generating Lund Jet Plane

### `lund-jet-generation.jl`

This script:
- Loads events from a sample data file
- Clusters final state particles using the anti-kt algorithm
- Sorts jets by transverse momentum squared
- Applies `generate_lund_projection` on each jet to calculate all the lund plane projections
- Prints the number of lund plane projections for each jet

#### How to Run

```julia
julia --project lund-jet-generation.jl
```

#### Customization

- **Event selection**: Change `event_no = 1` to select a different event.
- **Input data**: Replace the default path in `input_file` with your custom dataset.

---

### `lund-plane-visualisation.jl`

This script:
- Runs over all events in the dataset
- Extracts Lund-plane variables for each jet (`kt`, `ΔR`)
- Computes the average Lund plane using `generate_average_lund_image`
- Visualizes it using `CairoMakie`
- Saves the heatmap as `lund-gen.png`

#### How to Run

```julia
julia --project lund-plane-visualisation.jl
```

#### Output

- A PNG image `lund-gen.png` will be saved.

#### Customization

- **Input data**: Replace the default path in `input_file` with your custom dataset.
- **Jet algorithm parameters**: Modify `R = 1.0`, `ptmin = 10.0`, or the clustering algorithm.
- **Binning**: Adjust `bins` in `generate_average_lund_image` to control resolution.
- **Axis ranges**: Tune `xrange` and `yrange` to zoom in/out of particular regions.
- **Visuals**: Customize colormap or labels in the `heatmap!` function.

---
