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
struct ZeroRateCurve{num1 <: Number, num2 <: Number} <: AbstractZeroRateCurve
    r::CurveType{num1, num2}
    function ZeroRateCurve(r_::Array{num1}, T::num2) where {num1 <: Number, num2 <: Number}
        new{num2, num1}(Curve(r_, T))
    end

    function ZeroRateCurve(r_::CurveType{num1, num2}) where {num1 <: Number, num2 <: Number}
        new{num2, num1}(r_)
    end
end

#Common Constructors
function ZeroRate(r_::Array{num1}, T::num2) where {num1 <: Number, num2 <: Number}
    return ZeroRateCurve(r_, T)
end
"""
    ImpliedZeroRate(r::Vector, T::Vector)
    ImpliedZeroRate(r::Vector, T)

[`ZeroRateCurve`](@ref) reproducing the zero rates `r` quoted at maturities `T` (or at a uniform
grid up to `T`), i.e. whose instantaneous rate integrates to ``r_i T_i`` up to each ``T_i``.
See [`FinancialMonteCarlo.ImpliedCurve`](@ref).
"""
function ImpliedZeroRate(r_::Array{num1}, T::num2) where {num1 <: Number, num2 <: Number}
    return ZeroRateCurve(ImpliedCurve(r_, T))
end
function ImpliedZeroRate(r_::Array{num1}, T::Array{num2}) where {num1 <: Number, num2 <: Number}
    return ZeroRateCurve(ImpliedCurve(r_, T))
end

export ZeroRateCurve;
export ZeroRate;
