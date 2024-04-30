using EnumX

# Valid strategy enum (this is a scoped enum)
@enumx T=Strategy JetRecoStrategy Best N2Plain N2Tiled

const AllJetRecoStrategies = [ String(Symbol(x)) for x in instances(JetRecoStrategy.Strategy) ]

# Algorithm emun
@enumx T=Algorithm JetAlgorithm AntiKt Cambridge Kt EEKt Durham

const AllJetRecoAlgorithms = [ String(Symbol(x)) for x in instances(JetAlgorithm.Algorithm) ]

# Map from string to an enum value (used for CLI parsing)
Base.tryparse(E::Type{<:Enum}, str::String) =
        let insts = instances(E) ,
            p = findfirst(==(Symbol(str)) ∘ Symbol, insts) ;
            p !== nothing ? insts[p] : nothing
        end
