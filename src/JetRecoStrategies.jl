using MLStyle

# Valid strategy enum
@enum JetRecoStrategy Best N2Plain N2Tiled

# Map from string to strategy
Base.tryparse(E::Type{<:Enum}, str::String) =
        let insts = instances(E) ,
            p = findfirst(==(Symbol(str)) ∘ Symbol, insts) ;
            p !== nothing ? insts[p] : nothing
        end

const AllJetRecoStrategies = [ String(Symbol(x)) for x in instances(JetRecoStrategy) ]

# We will use the MLStyle package for a "switch", so needs a few tweaks for Julia enums
# https://thautwarm.github.io/MLStyle.jl/latest/syntax/pattern.html#support-pattern-matching-for-julia-enums
MLStyle.is_enum(::JetRecoStrategy) = true
MLStyle.enum_matcher(e::JetRecoStrategy, expr) = :($e === $expr)
