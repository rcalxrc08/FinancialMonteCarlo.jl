"""
Struct for Geometric Brownian Motion

		gbmProcess=GeometricBrownianMotion(σ::num1,μ::num2) where {num1,num2 <: Number}
	
Where:\n
		σ	=	volatility of the process.
		μ	=	drift of the process.
		x0	=	initial value.
"""
mutable struct GeometricBrownianMotion{num <: Number, num1 <: Number, num2 <: Number, numtype <: Number} <: ItoProcess{numtype}
    σ::num
    μ::num1
    x0::num2
    function GeometricBrownianMotion(σ::num, μ::num1, x0::num2) where {num <: Number, num1 <: Number, num2 <: Number}
        if σ <= 0
            error("Volatility must be positive")
        else
            zero_typed = zero(num) + zero(num1)
            return new{num, num1, num2, typeof(zero_typed)}(σ, μ, x0)
        end
    end
end

export GeometricBrownianMotion;

function simulate!(X, mcProcess::GeometricBrownianMotion, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    @assert T > 0
    σ_gbm = mcProcess.σ
    mu_gbm = mcProcess.μ
    μ_bm = mu_gbm - σ_gbm^2 / 2
    simulate!(X, BrownianMotion(σ_gbm, μ_bm), mcBaseData, T)
    S0 = mcProcess.x0
    # @. X=S0*exp(X);
    broadcast!(x -> S0 * exp(x), X, X)

    nothing
end
