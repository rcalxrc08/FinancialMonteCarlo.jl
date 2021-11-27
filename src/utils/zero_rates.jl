###########Zero Rate Curve Types
#Abstract
abstract type AbstractZeroRateCurve end

#Scalar
struct ZeroRate{T2 <: Number} <: AbstractZeroRateCurve
    r::T2
    function ZeroRate(r::T2) where {T2 <: Number}
        return new{T2}(r)
    end
end

#Curve
struct ZeroRateCurve{num1 <: Number, num2 <: Number} <: AbstractZeroRateCurve
    r::Curve{num1, num2}
    function ZeroRateCurve(r_::Array{num1}, T::num2) where {num1 <: Number, num2 <: Number}
        new{num2, num1}(Curve(r_, T))
    end

    function ZeroRateCurve(r_::Curve{num1, num2}) where {num1 <: Number, num2 <: Number}
        new{num2, num1}(r_)
    end
end

#Common Constructors
function ZeroRate(r_::Array{num1}, T::num2) where {num1 <: Number, num2 <: Number}
    return ZeroRateCurve(r_, T)
end
function ImpliedZeroRate(r_::Array{num1}, T::num2) where {num1 <: Number, num2 <: Number}
    return ZeroRateCurve(ImpliedCurve(r_, T))
end
function ImpliedZeroRate(r_::Array{num1}, T::Array{num2}) where {num1 <: Number, num2 <: Number}
    return ZeroRateCurve(ImpliedCurve(r_, T))
end

export ZeroRateCurve;
export ZeroRate;
