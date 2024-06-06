# EnumX provides scoped enums, which are nicer
using EnumX

# Valid strategy enum
@enumx T = Strategy RecoStrategy Best N2Plain N2Tiled
const AllJetRecoStrategies = [String(Symbol(x)) for x in instances(RecoStrategy.Strategy)]

# Algorithm emun
@enumx T = Algorithm JetAlgorithm AntiKt CA Kt EEKt Durham
const AllJetRecoAlgorithms = [String(Symbol(x)) for x in instances(JetAlgorithm.Algorithm)]

# Map from algorithms to power values used 
power2algorithm = Dict(-1 => JetAlgorithm.AntiKt,
    0 => JetAlgorithm.CA,
    1 => JetAlgorithm.Kt)
algorithm2power = Dict(JetAlgorithm.AntiKt => -1,
    JetAlgorithm.CA => 0,
    JetAlgorithm.Kt => 1)

# Map from string to an enum value (used for CLI parsing)
Base.tryparse(E::Type{<:Enum}, str::String) =
    let insts = instances(E),
        p = findfirst(==(Symbol(str)) ∘ Symbol, insts)

        p !== nothing ? insts[p] : nothing
    end
