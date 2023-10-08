# Abstract type for GeometricBrownianMotion, an Ito process is always a Levy
abstract type AbstractGeometricBrownianMotion{T <: Number} <: AbstractMonteCarloEngine{T} end
"""
Struct for Geometric Brownian Motion

		gbmProcess=GeometricBrownianMotion(σ::num1,μ::num2) where {num1,num2 <: Number}
	
Where:

		σ	=	volatility of the process.
		μ	=	drift of the process.
		x0	=	initial value.
"""
struct GeometricBrownianMotion{num <: Number, num1 <: Number, num2 <: Number, numtype <: Number} <: AbstractGeometricBrownianMotion{numtype}
    σ::num
    μ::num1
    x0::num2
    function GeometricBrownianMotion(σ::num, μ::num1, x0::num2) where {num <: Number, num1 <: Number, num2 <: Number}
        ChainRulesCore.@ignore_derivatives @assert σ > 0 "Volatility must be positive"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num) + zero(num1)
        return new{num, num1, num2, typeof(zero_typed)}(σ, μ, x0)
    end
end

export GeometricBrownianMotion;

function simulate!(X, mcProcess::AbstractGeometricBrownianMotion, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    ChainRulesCore.@ignore_derivatives @assert T > 0
    σ = mcProcess.σ
    μ = mcProcess.μ - σ^2 / 2
    simulate!(X, BrownianMotion(σ, μ), mcBaseData, T)
    S0 = mcProcess.x0
    @. X = S0 * exp(X)

    nothing
end

# function simulate_path!(X, mcProcess::AbstractGeometricBrownianMotion, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
# ChainRulesCore.@ignore_derivatives @assert T > 0
#     σ = mcProcess.σ
#     μ = mcProcess.μ - σ^2 / 2
#     simulate_path!(X, BrownianMotion(σ, μ), mcBaseData, T)
#     S0 = mcProcess.x0
#     @. X = S0 * exp(X)

#     nothing
# end