"""Interface for composite types that includes fields px, py, py, and E
that represents the components of a four-momentum vector."""
abstract type FourMomentum end

"""
Simple jet structure without cluster history
"""
struct PlainJet <: FourMomentum
    px::Float64
    py::Float64
    pz::Float64
    E::Float64
end
