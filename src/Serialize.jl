"""
`savejets(filename, jets; format="px py pz E")`

Saves the given `jets` into a file with the given `filename`. Each line contains information about a single jet and is formatted according to the `format` string which defaults to `"px py pz E"` but can also contain other values in any order: `"pt2"` for pt^2, `"phi"` for azimuth, `"rapidity"` for rapidity. It is strongly NOT recommended to put something other than values and (possibly custom) separators in the `format` string.
"""
function savejets(filename, jets; format="px py pz E")
    symbols = Dict(
        "E" => JetReconstruction.energy,
        "energy" => JetReconstruction.energy,
        "px" => JetReconstruction.px,
        "py" => JetReconstruction.py,
        "pz" => JetReconstruction.pz,
        "pt2" => JetReconstruction.pt2,
        "phi" => JetReconstruction.phi,
        "rapidity" => JetReconstruction.rapidity
    )
    for pair in symbols
        if !occursin(pair[1], format)
            pop!(symbols, pair[1])
        end
    end

    open(filename, "w") do file
        write(file, "# this file contains jet data.\n# each line that does not start with '#' contains the information about a jet in the following format:\n# "*"$(format)"*"\n# where E is energy, px is momentum along x, py is momentum along y, pz is momentum along z, pt2 is pt^2, phi is azimuth, and rapidity is rapidity\n")
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
    loadjets!(filename, jets; splitby=isspace, constructor=(x,y,z,E)->LorentzVector(E,px,py,pz), dtype=Float64) -> jets

Loads the `jets` from a file. Ignores lines that start with `'#'`. Each line gets processed in the following way: the line is split using `split(line, splitby)` or simply `split(line)` by default. Every value in this line is then converted to the `dtype` (which is `Float64` by default). These values are then used as arguments for the `constructor` function which should produce individual jets. By default, the `constructor` constructs Lorentz vectors.

Everything that was already in `jets` is not affected as we only use `push!` on it.
```julia
# example that loads jets from two files into one array
jets = []
loadjets!("myjets1.dat", jets)
loadjets!("myjets2.dat", jets)
```
"""
function loadjets!(filename, jets; splitby=isspace, constructor=(x,y,z,E)->LorentzVector(E,px,py,pz), dtype=Float64)
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
    loadjets(filename; splitby=isspace, constructor=(px,py,pz,E)->LorentzVector(E,px,py,pz), VT=LorentzVector) -> jets

Loads the `jets` from a file, where each element of `jets` is of type `VT`. Ignores lines that start with `'#'`. Each line gets processed in the following way: the line is split using `split(line, splitby)` or simply `split(line)` by default. These values are then used as arguments for the `constructor` function which should produce individual jets of the `VT` type. By default, the `constructor` constructs Lorentz vectors.

```julia
# example
jets = loadjets("myjets1.dat")
```
"""
function loadjets(filename; splitby=isspace, constructor=(px,py,pz,E)->LorentzVector(E,px,py,pz), VT=LorentzVector)
    loadjets!(filename, VT[], splitby=splitby, constructor=constructor)
end
