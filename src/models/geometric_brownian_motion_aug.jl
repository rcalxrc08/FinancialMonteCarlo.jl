"""
Struct for Geometric Brownian Motion

		gbmProcess=GeometricBrownianMotionVec(σ::num1,μ::num2) where {num1,num2 <: Number}
	
Where:

		σ	=	volatility of the process.
		μ	=	drift of the process.
"""
mutable struct GeometricBrownianMotionVec{num <: Number, num1 <: Number, num4 <: Number, num2 <: Number, numtype <: Number} <: AbstractGeometricBrownianMotion{numtype}
    σ::num
    μ::CurveType{num1, num4}
    x0::num2
    function GeometricBrownianMotionVec(σ::num, μ::CurveType{num1, num4}, x0::num2) where {num <: Number, num1 <: Number, num4 <: Number, num2 <: Number}
        @assert σ > 0 "Volatility must be positive"
        zero_typed = zero(num) + zero(num1) + zero(num4) + zero(num2)
        return new{num, num1, num4, num2, typeof(zero_typed)}(σ, μ, x0)
    end
end

function GeometricBrownianMotion(σ::num, μ::CurveType{num1, num4}, x0::num2) where {num <: Number, num1 <: Number, num4 <: Number, num2 <: Number}
    return GeometricBrownianMotionVec(σ, μ, x0)
end

export GeometricBrownianMotionVec;