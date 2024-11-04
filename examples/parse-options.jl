"""
Add `parse_item` code for interpreting `JetAlgorithm.Algorithm` and `RecoStrategy.Strategy` types
from the command line.
"""

function ArgParse.parse_item(::Type{RecoStrategy.Strategy}, x::AbstractString)
    s = tryparse(RecoStrategy.Strategy, x)
    if s === nothing
        throw(ErrorException("Invalid value for strategy: $(x)"))
    end
    s
end

function ArgParse.parse_item(::Type{JetAlgorithm.Algorithm}, x::AbstractString)
    s = tryparse(JetAlgorithm.Algorithm, x)
    if s === nothing
        throw(ErrorException("Invalid value for algorithm: $(x)"))
    end
    s
end
