"""
Struct for LogNormalMixture

		lnmModel=LogNormalMixture(η::Array{num1},λ::Array{num2}) where {num1,num2<: Number}
	
Where:\n
		η  =	Array of volatilities.
		λ  = 	Array of weights.
"""
mutable struct LogNormalMixture{num <: Number, num2 <: Number, abstrUnderlying <: AbstractUnderlying, numtype <: Number} <: ItoProcess{numtype}
    η::Array{num, 1}
    λ::Array{num2, 1}
    underlying::abstrUnderlying
    function LogNormalMixture(η::Array{num, 1}, λ::Array{num2, 1}, underlying::abstrUnderlying) where {num <: Number, num2 <: Number, abstrUnderlying <: AbstractUnderlying}
        @assert minimum(η) > 0.0 "Volatilities must be positive"
        @assert minimum(λ) > 0.0 "weights must be positive"
        @assert sum(λ) <= 1.0 "λs must be weights"
        @assert length(λ) == length(η) - 1 "Check vector lengths"
        zero_typed = zero(num) + zero(num2)
        return new{num, num2, abstrUnderlying, typeof(zero_typed)}(η, λ, underlying)
    end
end

export LogNormalMixture;

function simulate!(S, mcProcess::LogNormalMixture, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    @assert T > 0.0
    r = rfCurve.r
    S0 = mcProcess.underlying.S0
    d = dividend(mcProcess)
    η_gbm = mcProcess.η
    λ_gmb = mcProcess.λ
    push!(λ_gmb, 1.0 - sum(mcProcess.λ))
    mu_gbm = r - d
    simulate!(S, GeometricBrownianMotion(η_gbm[1], mu_gbm, S0), mcBaseData, T)
    S .*= λ_gmb[1]
    tmp_s = similar(S)
    for (η_gbm_, λ_gmb_) in zip(η_gbm[2:end], λ_gmb[2:end])
        simulate!(tmp_s, GeometricBrownianMotion(η_gbm_, mu_gbm, S0), mcBaseData, T)
        S .+= λ_gmb_ * tmp_s
    end
    return nothing
end
