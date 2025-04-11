"""
    abstract type FourMomentum end

Interface for composite types that includes fields px, py, py, and E
that represents the components of a four-momentum vector. All concrete
jets that are used in the package are subtypes of this type.
"""
abstract type FourMomentum end


"""
    Jet <: FourMomentum

An abstract type representing a jet, rather than just a plain FourMomentum.

# Details
This allows for more efficient extraction of certain properties of the jet.
"""
abstract type Jet <: FourMomentum end


"""
    struct PlainJet <: FourMomentum

Simple jet structure without cluster history, used as an intermediate
type during jet merging.
"""
struct PlainJet <: FourMomentum
    px::Float64
    py::Float64
    pz::Float64
    E::Float64
end

# Define here all common functions that can be used for all jet types <: FourMomentum
import Base.+;
"""
    +(jet1::FourMomentum, jet2::FourMomentum)

Adds two four-momentum vectors together, returning a new `PlainJet` jet.
"""
function +(jet1::FourMomentum, jet2::FourMomentum)
    PlainJet(jet1.px + jet2.px, jet1.py + jet2.py, jet1.pz + jet2.pz, jet1.E + jet2.E)
end

px(j::FourMomentum) = j.px
py(j::FourMomentum) = j.py
pz(j::FourMomentum) = j.pz
energy(j::FourMomentum) = j.E
p2(j::FourMomentum) = j.px^2 + j.py^2 + j.pz^2
pt2(j::FourMomentum) = j.px^2 + j.py^2
pt(j::FourMomentum) = sqrt(pt2(j))

const kt2 = pt2
const E = energy

# Define here all common functions that can be used for all jet types <: Jet
cluster_hist_index(j::Jet) = j._cluster_hist_index
