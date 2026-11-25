###########Underlying Types
#Abstract
abstract type AbstractUnderlying end

"""
    Underlying(S0, d = zero(S0))

Spot price `S0` and continuous dividend yield `d` of an asset, passed next to the model to every entry point
(`pricer(model, underlying, zeroCurve, method, payoff)`): the models hold only their parameters (#72).
`d` is either a number or a time dependent curve (see [`FinancialMonteCarlo.Curve`](@ref) and
[`FinancialMonteCarlo.ImpliedCurve`](@ref)): one type, whose methods dispatch on the type of `d` where they differ.
`S0` must be non negative.

# Examples
```julia
Underlying(100.0)
Underlying(100.0, 0.01)
Underlying(100.0, FinancialMonteCarlo.Curve([0.01, 0.02], 1.0))
```
"""
struct Underlying{S <: Number, D <: Union{Number, AbstractCurve}} <: AbstractUnderlying
    S0::S
    d::D
    function Underlying(S0::S, d::D = zero(S0)) where {S <: Number, D <: Union{Number, AbstractCurve}}
        check_spot(S0)
        return new{S, D}(S0, d)
    end
end

# The continuous dividend yield of an underlying, a number or a curve (#72: read from the underlying passed next to the
# model, not from the model)
dividend(und::AbstractUnderlying) = und.d

# The underlyings of a multivariate model: one per component, in the order of its models
const UnderlyingTuple = Tuple{AbstractUnderlying, Vararg{AbstractUnderlying}}
# What an entry point takes next to the model
const UnderlyingArg = Union{AbstractUnderlying, UnderlyingTuple}

export Underlying;
