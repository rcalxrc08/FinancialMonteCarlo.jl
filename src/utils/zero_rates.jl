###########Zero Rate Curve Types
#Abstract
abstract type AbstractZeroRateCurve end

#Scalar
"""
    ZeroRate(r)
    ZeroRate(r::Vector, T)

Zero rate curve used for discounting and for the risk neutral drift.
With a number, the rate is constant; with a vector, `r` holds the instantaneous rates on a uniform
grid from 0 to `T` and a [`ZeroRateCurve`](@ref) is returned (see [`FinancialMonteCarlo.Curve`](@ref)).
To build the curve from quoted zero rates use [`FinancialMonteCarlo.ImpliedZeroRate`](@ref).
"""
struct ZeroRate{T2 <: Number} <: AbstractZeroRateCurve
    r::T2
    function ZeroRate(r::T2) where {T2 <: Number}
        return new{T2}(r)
    end
end

#Curve
"""
    ZeroRateCurve(r::Vector, T)
    ZeroRateCurve(curve)

Time dependent zero rate curve, holding a curve of instantaneous rates (see
[`FinancialMonteCarlo.Curve`](@ref)). Usually built through [`ZeroRate`](@ref) or
[`FinancialMonteCarlo.ImpliedZeroRate`](@ref).
"""
struct ZeroRateCurve{C <: Curve} <: AbstractZeroRateCurve
    r::C
    ZeroRateCurve(r_::C) where {C <: Curve} = new{C}(r_)
end
ZeroRateCurve(r_::Union{AbstractVector{<:Number}, Tuple{Vararg{Number}}}, T::Number) = ZeroRateCurve(Curve(r_, T))

#Common Constructors
function ZeroRate(r_::Union{AbstractVector{<:Number}, Tuple{Vararg{Number}}}, T::Number)
    return ZeroRateCurve(r_, T)
end
"""
    ImpliedZeroRate(r::Vector, T::Vector)
    ImpliedZeroRate(r::Vector, T)

[`ZeroRateCurve`](@ref) reproducing the zero rates `r` quoted at maturities `T` (or at a uniform
grid up to `T`), i.e. whose instantaneous rate integrates to ``r_i T_i`` up to each ``T_i``.
See [`FinancialMonteCarlo.ImpliedCurve`](@ref).
"""
ImpliedZeroRate(r_, T) = ZeroRateCurve(ImpliedCurve(r_, T))

export ZeroRateCurve;
export ZeroRate;
