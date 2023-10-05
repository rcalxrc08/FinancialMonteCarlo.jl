"""
Struct for VG Process

		vgProcess=VarianceGammaProcess(σ::num1,θ::num2,κ::num3) where {num1,num2,num3<: Number}
	
Where:

		σ =	volatility of the process.
		θ = variance of volatility.
		κ =	skewness of volatility.
"""
struct VarianceGammaProcess{num <: Number, num1 <: Number, num2 <: Number, abstrUnderlying <: AbstractUnderlying, numtype <: Number} <: InfiniteActivityProcess{numtype}
    σ::num
    θ::num1
    κ::num2
    underlying::abstrUnderlying
    function VarianceGammaProcess(σ::num, θ::num1, κ::num2, underlying::abstrUnderlying) where {num <: Number, num1 <: Number, num2 <: Number, abstrUnderlying <: AbstractUnderlying}
        @assert σ > 0 "volatility must be positive"
        @assert κ > 0 "κappa must be positive"
        @assert 1 - σ * σ * κ / 2 - θ * κ >= 0 "Parameters with unfeasible values"
        zero_typed = zero(num) + zero(num1) + zero(num2)
        return new{num, num1, num2, abstrUnderlying, typeof(zero_typed)}(σ, θ, κ, underlying)
    end
end

export VarianceGammaProcess;

function simulate!(X, mcProcess::VarianceGammaProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    r = rfCurve.r
    S0 = mcProcess.underlying.S0
    d = dividend(mcProcess)
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    θ = mcProcess.θ
    κ = mcProcess.κ
    @assert T > 0
    dt = T / Nstep
    ψ = -1 / κ * log(1 - σ^2 * κ / 2 - θ * κ)
    drift = r - d - ψ

    gammaRandomVariable = Gamma(dt / κ, κ)

    simulate!(X, SubordinatedBrownianMotion(σ, drift, θ, gammaRandomVariable), mcBaseData, T)

    @. X = S0 * exp(X)

    nothing
end
