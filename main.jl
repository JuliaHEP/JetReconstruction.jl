# A temporary file to run some quick tests

using StaticArrays

include("src/JetReconstruction.jl")
using .JetReconstruction
# using JetReconstruction

## Unrealistic tests
anti_kt([
    SVector(π, 0, 0, 1)
])

anti_kt([
    SVector(π, 0, 0, 1),
    SVector(π, 0, 2, 0),
    SVector(1, 1, 0, 0),
    SVector(1, 0, 1, 0),
])

## Developer convenience test (running the algo on a custom data structure)
import .JetReconstruction # no need to import Particle

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
    return CylVector(pt,eta,phi,mass)
end

@time anti_kt([
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
@time anti_kt([
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
