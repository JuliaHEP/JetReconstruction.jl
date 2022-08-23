"""
Defines `savejets`, `loadjets`, `loadjets!`.
"""
module Serialize

import ..JetReconstruction

export savejets, loadjets!, loadjets

"""
`savejets(filename, jets; format="px py pz E")`

Saves the given `jets` into a file with the given `filename`. Each line contains information about a single jet and is formatted according to the `format` string which defaults to `"E px py pz"` but can also contain other values in any order: `"pt"` or `"kt"` for transverse momentum, `"phi"` for azimuth, `"eta"` for pseudorapidity, `"m"` for mass. It is strongly NOT recommended to put something other than values and (possibly custom) separators in the `format` string.
"""
function savejets(filename, jets; format="px py pz E")
    symbols = Dict(
        "E" => JetReconstruction.energy,
        "px" => JetReconstruction.px,
        "pt" => JetReconstruction.pt,
        "kt" => JetReconstruction.pt,
        "py" => JetReconstruction.py,
        "pz" => JetReconstruction.pz,
        "phi" => JetReconstruction.phi,
        "m" => JetReconstruction.mass,
        "eta" => JetReconstruction.eta
    )
    for pair in symbols
        if !occursin(pair[1], format)
            pop!(symbols, pair[1])
        end
    end

    open(filename, "w") do file
        write(file, "# this file contains jet data.\n# each line that does not start with '#' contains the information about a jet in the following format:\n# "*"$(format)"*"\n# where E is energy, px is momentum along x, py is momentum along y, pz is momentum along z, pt or kt is transverse momentum, phi is azimuth, eta is pseudorapidity, m is mass\n")
        for j in jets
            line = format
            for pair in symbols
                line = replace(line, pair[1] => "$(pair[2](j))")
            end
            write(file, line*'\n')
        end
        write(file, "#END")
    end
    nothing
end

"""
`loadjets!(filename, jets; splitby=isspace, constructor=(px,py,pz,E)->[px,py,pz,E], dtype=Float64) -> jets`

Loads the `jets` from a file. Ignores lines that start with `'#'`. Each line gets processed in the following way: the line is split using `split(line, splitby)` or simply `split(line)` by default. Every value in this line is then converted to the `dtype` (which is `Float64` by default). These values are then used as arguments for the `constructor` function which should produce individual jets. By default, the `constructor` constructs vectors of the form `[px,py,pz,E]`.

Everything that was already in `jets` is not affected as we only use `push!` on it.
```julia
# example that loads jets from two files into one array
jets = []
loadjets!("myjets1.dat", jets)
loadjets!("myjets2.dat", jets)
```
"""
function loadjets!(filename, jets; splitby=isspace, constructor=(x,y,z,E)->Float64[x,y,z,E], dtype=Float64)
    open(filename, "r") do file
        for line in eachline(file)
            if line[1] != '#'
                jet = constructor(
                    (parse(dtype, x) for x in split(line, splitby) if x != "")...
                )
                push!(jets, jet)
            end
        end
    end

    jets
end

"""
`loadjets(filename; splitby=isspace, constructor=(px,py,pz,E)->Float64[px,py,pz,E], dtype=Float64) -> jets`

Loads the `jets` from a file. Ignores lines that start with `'#'`. Each line gets processed in the following way: the line is split using `split(line, splitby)` or simply `split(line)` by default. Every value in this line is then converted to the `dtype` (which is `Float64` by default). These values are then used as arguments for the `constructor` function which should produce individual jets. By default, the `constructor` constructs vectors of the form `[px,py,pz,E]`.

```julia
# example
jets = loadjets("myjets1.dat")
```
"""
function loadjets(filename; splitby=isspace, constructor=(px,py,pz,E)->Float64[px,py,pz,E], dtype=Float64)
    loadjets!(filename, [], splitby=splitby, constructor=constructor)
end

end
