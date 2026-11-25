###########Zero Rate Curve Types
#Abstract
abstract type AbstractZeroRateCurve end

"""
    ZeroRate(r)

Zero rate curve used for discounting and for the risk neutral drift: `r` is a constant rate, or a curve of
instantaneous rates (a [`FinancialMonteCarlo.Curve`](@ref), or a [`FinancialMonteCarlo.ImpliedCurve`](@ref) to
reproduce quoted zero rates). One type, whose methods dispatch on the type of `r` where they differ.

# Examples
```julia
ZeroRate(0.02)
ZeroRate(FinancialMonteCarlo.Curve([0.01, 0.02], 1.0))          # instantaneous rates on a uniform grid up to 1
ZeroRate(FinancialMonteCarlo.ImpliedCurve([0.01, 0.02], [0.5, 1.0]))   # zero rates quoted at 0.5 and 1
```
"""
struct ZeroRate{R <: Union{Number, AbstractCurve}} <: AbstractZeroRateCurve
    r::R
    ZeroRate(r::R) where {R <: Union{Number, AbstractCurve}} = new{R}(r)
end

export ZeroRate;
