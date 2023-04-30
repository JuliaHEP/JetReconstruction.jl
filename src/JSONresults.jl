using JSON3
using StructTypes

"Simple structures to capture the state of a final jet output
from the algorithm"
struct FinalJet
    rap::Float64
    phi::Float64
    pt::Float64
end

struct FinalJets
    jetid::Int64
    jets::Vector{FinalJet}
end

# Register types to be able to dump to JSON
StructTypes.StructType(::Type{FinalJet}) = StructTypes.Struct()
StructTypes.StructType(::Type{FinalJets}) = StructTypes.Struct()
