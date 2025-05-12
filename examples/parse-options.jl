"""
Add `parse_item` code for interpreting `JetAlgorithm.Algorithm` and
`RecoStrategy.Strategy` types from the command line.
"""

function do_enum_parse(E::Type, x::AbstractString)
    insts = instances(E)
    p = findfirst(==(Symbol(x)) ∘ Symbol, insts)
    p !== nothing ? insts[p] : nothing
end

function ArgParse.parse_item(E::Type{JetAlgorithm.Algorithm}, x::AbstractString)
    p = do_enum_parse(E, x)
    if p === nothing
        throw(ErrorException("Invalid value for algorithm: $(x)"))
    end
    p
end

function ArgParse.parse_item(E::Type{RecoStrategy.Strategy}, x::AbstractString)
    p = do_enum_parse(E, x)
    if p === nothing
        throw(ErrorException("Invalid value for strategy: $(x)"))
    end
    p
end

function ArgParse.parse_item(E::Type{RecombinationScheme.Recombine}, x::AbstractString)
    p = do_enum_parse(E, x)
    if p === nothing
        throw(ErrorException("Invalid value for recombination scheme: $(x)"))
    end
    p
end
