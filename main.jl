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

using CairoMakie

## array constructor
makejet(pt, ϕ, y) = [pt*cos(ϕ), pt*sin(ϕ), pt*sinh(y), sqrt((pt*cos(ϕ))^2 + (pt*sin(ϕ))^2 + (pt*sinh(y))^2)]

## choose a data sample and save it
smalldata = [
    makejet(25, 1.2, 0),
    makejet(0.1, 1, 0),
    makejet(20, 0, 0)
]
smalldata = [
    makejet(25, 1.2, 0),
    makejet(0.1, 0.6, 0),
    makejet(20, 0, 0)
]
smalldata = [
    makejet(5, 0.8, 0),
    makejet(2, 0, 0.8),
]
smalldata = [
    makejet(30, -0.5, -2.5),
    makejet(30, -0.5, -3.5)
]
smalldata = [
    makejet(30, 0.5, -2.5),
    makejet(30, 0.5, -3.5),
    makejet(30, -0.5, -2.5),
    makejet(30, -0.5, -3.5),
    makejet(0.1, 0, -3),
]
smalldata = [
    [0.2839991726, -0.4668202293, -49.6978846131, 49.709744138],
    [0.1434555273, -0.0312880942, -6.9382351613, 6.9411919273],
]
smalldata = [
    makejet(0.8, 0.8, 0.1),
    makejet(0.1, 0.1, 0.8),
]

smalldata = [
    [-0.8807412236, -1.2331262152, -157.431315651, 157.4393822839],
    [-0.0051712611, 0.23815508, -9.7396045662, 9.7435168946], # 2
    [0.0362280943, 0.2694752057, -6.9243427525, 6.9310844534], # 3
    [-0.2206628664, -0.1438198985, -0.6838608666, 0.7460038429],
    [1.2716787521, 1.0422298083, -6.1740167274, 6.3907254797], # 5
    [-0.5695590845, -0.3627761836, -58.5430479911, 58.5544811606],
    [0.2839991726, -0.4668202293, -49.6978846131, 49.709744138],
    [0.6510530003, 1.3970949413, -62.7226079598, 62.7485783532], # 8
    [0.1434555273, -0.0312880942, -6.9382351613, 6.9411919273],
    [0.4931562547, 2.1627817414, -14.8865871635, 15.0516040711], # 10
    [0.2396813608, -0.0786236784, -1.9340954697, 1.9554625817],
    [0.3355486441, 0.0516402769, -0.8346540063, 0.9118040941],
    [-0.7853865645, -0.7810520475, -1.5367790662, 1.8994852039],
]
smalldata = [
    [1.2716787521, 1.0422298083, -6.1740167274, 6.3907254797],
    [0.4931562547, 2.1627817414, -14.8865871635, 15.0516040711],
    [0.3355486441, 0.0516402769, -0.8346540063, 0.9118040941],
]

savejets("small.dat", smalldata, format="px py pz E")
## run
smalljets, smallind = anti_kt_algo(smalldata, R=1); img = jetsplot(smalldata, smallind, Module=CairoMakie)

smalljets, smallind = anti_kt_algo_alt(smalldata, R=1); img = jetsplot(smalldata, smallind, Module=CairoMakie)

smalljets, smallind = kt_algo(smalldata, R=1); img = jetsplot(smalldata, smallind, Module=CairoMakie)
