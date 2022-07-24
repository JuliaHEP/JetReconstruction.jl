module Serialize

import ..JetReconstruction

export savejets, loadjets

"""
`savejets(filename, jets; format="E px py pz")`

Saves the given `jets` into a file with the given `filename`. Each line contains information about a single jet and is formatted according to the `format` string which defaults to `"E px py pz"` but can also contain other values in any order: `"pt"` or `"kt"` for transverse momentum, `"phi"` for azimuth, `"eta"` for pseudorapidity, `"m"` for mass.
"""
function savejets(filename, jets; format="E px py pz")
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

end
