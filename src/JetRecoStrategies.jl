using EnumX

# Valid strategy enum (this is a scoped enum)
@enumx JetRecoStrategy Best N2Plain N2Tiled

# Map from string to an enum value (used for CLI parsing)
Base.tryparse(E::Type{<:Enum}, str::String) =
        let insts = instances(E) ,
            p = findfirst(==(Symbol(str)) ∘ Symbol, insts) ;
            p !== nothing ? insts[p] : nothing
        end

const AllJetRecoStrategies = [ String(Symbol(x)) for x in instances(JetRecoStrategy.T) ]

