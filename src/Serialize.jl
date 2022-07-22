"""
`savejets(filename, jets; format="E px py pz")`

Saves the given `jets` into a file with the given `filename`. Each line contains information about a single jet and is formatted according to the `format` string which defaults to `"E px py pz"` but can also contain other values in any order: `"pt"` or `"kt"` for transverse momentum, `"phi"` for azimuth, `"eta"` for pseudorapidity, `"m"` for mass.
"""
function savejets(filename, jets; format="E px py pz")
    open(filename, "w") do file
        write(file, "# this file contains jet data.\n# each line that does not start with '#' contains the information about a jet in the following format:\n# "*"$(format)"*"\n# where E is energy, px is momentum along x, py is momentum along y, pz is momentum along z, pt or kt is transverse momentum, phi is azimuth, eta is pseudorapidity, m is mass\n")
        for j in jets
            line = format
            try
                line = replace(line, "E" => "$(Main.energy(j))")
            catch; end
            try
            jpt = "$(Main.pt(j))"
            line = replace(line, "pt" => jpt)
            line = replace(line, "kt" => jpt)
            catch; end
            try
            line = replace(line, "phi" => "$(Main.phi(j))")
            catch; end
            try
            line = replace(line, "eta" => "$(Main.eta(j))")
            catch; end
            try
            line = replace(line, "m" => "$(Main.mass(j))")
            catch; end
            try
            line = replace(line, "px" => "$(Main.px(j))")
            catch; end
            try
            line = replace(line, "py" => "$(Main.py(j))")
            catch; end
            try
            line = replace(line, "pz" => "$(Main.pz(j))")
            catch; end
            write(file, line*'\n')
        end
        write(file, "#END")
    end
    nothing
end
