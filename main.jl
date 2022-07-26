# A temporary file to run some quick tests

# choose one:
# this
#=
include("src/JetReconstruction.jl")
using .JetReconstruction
=#
# or this
using Revise; import Pkg
Pkg.activate(".")
using JetReconstruction

## Realistic test
using StaticArrays
datalen = 10

data = [SVector{4, Float64}[] for _ in 1:datalen]
test_number = 1
for line in eachline("test/data/Pythia-PtMin1000-LHC-10ev.dat")
    if line == "#END"
        test_number += 1
    elseif line[1] != '#'
        px, py, pz, E = (parse(Float64, x) for x in split(line))
        vec = SVector(E, px, py, pz)
        push!(data[test_number], vec)
    end
end

for i in 1:datalen
    savejets("test/data/"*string(i)*".dat", data[i], format="px py pz E")
end

loadjets("test/data/1-fj-result.dat", constructor=(x,y,z,E)->[E,x,y,z])
anti_kt(data[1])[1]

precompile(anti_kt, (typeof(data[1]),))
precompile(anti_kt_alt, (typeof(data[1]),))

jetarrs = []
objectidxarrs = Vector{Vector{Int}}[]
for i in 1:1#datalen
    jetarr, components = anti_kt(data[i])
    softs = [j for j in 1:length(components) if (length(components[j]) == 1 || jetarr[j][1] < 2)]
    deleteat!(jetarr, softs)
    deleteat!(components, softs)
    push!(jetarrs, jetarr)
    push!(objectidxarrs, components)
end

## Save or load data
savejets("jetsavetest.dat", data[1])
somedata = loadjets("jetsavetest.dat", constructor=SVector)
somedata == data[1]

## Visualisation
index = 10
img = jetsplot(data[index], objectidxarrs[index])

#display(img) # for Juno Plots window
#PyPlot.show() # for terminal usage

## Developer convenience test (running the algo on a custom data structure)
import JetReconstruction # no need to import Particle

struct CylVector
    y::Float64
    ϕ::Float64
    pt::Float64
    mass::Float64
end
CylVector(y, ϕ, pt) = CylVector(y, ϕ, pt, 0)

JetReconstruction.eta(x::CylVector) = x.y
JetReconstruction.phi(x::CylVector) = x.ϕ
JetReconstruction.pt(x::CylVector) = x.pt
JetReconstruction.mass(x::CylVector) = x.mass
function Base.:+(x::CylVector, y::CylVector)
    px(v::CylVector) = v.pt * cos(v.ϕ)
    py(v::CylVector) = v.pt * sin(v.ϕ)
    pz(v::CylVector) = v.pt * sinh(v.y)
    energy(v::CylVector) = sqrt(px(v)^2 + py(v)^2 + pz(v)^2 + v.mass^2)
    m1, m2 = max(x.mass, 0), max(y.mass, 0)

    px1, px2 = px(x), px(y)
    py1, py2 = py(x), py(y)
    pz1, pz2 = pz(x), pz(y)
    e1 = sqrt(px1^2 + py1^2 + pz1^2 + m1^2)
    e2 = sqrt(px2^2 + py2^2 + pz2^2 + m2^2)

    sumpx = px1+px2
    sumpy = py1+py2
    sumpz = pz1+pz2

    ptsq = sumpx^2 + sumpy^2
    pt = sqrt(ptsq)
    eta = asinh(sumpz/pt)
    phi = atan(sumpy, sumpx)
    mass = sqrt(max(muladd(m1, m1, m2^2) + 2*e1*e2 - 2*(muladd(px1, px2, py1*py2) + pz1*pz2), 0))
    return CylVector(eta,phi,pt,mass)
end

cyljets, _ = @time anti_kt([
    CylVector(0, 0, 130),
    CylVector(0, 0.7, 120),
    CylVector(0, 0.7, 80),
    CylVector(0, 1.5, 90)
])
@time anti_kt([
    CylVector(0, 0, 130),
    CylVector(0, 0.7, 200),
    CylVector(0, 1.5, 90)
])
cyljets, _ = @time anti_kt([
    CylVector(0, 1, 130),
    CylVector(0, 0.7, 120),
    CylVector(0, 0.01, 80),
    CylVector(0, 0.01, 90),
    CylVector(0, 0.11, 81),
    CylVector(0, 0.02, 90),
    CylVector(0, 0.02, 80),
    CylVector(0, 0.03, 90),
    CylVector(0, 0.04, 83),
    CylVector(0, 0.05, 90)
])

savejets("jetsavetest2.dat", cyljets, format="eta phi kt m")
somedata2 = loadjets("jetsavetest2.dat", constructor=CylVector)
somedata2 == cyljets
