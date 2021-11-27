"""
Struct for ShiftedLogNormalMixture

		lnmModel=ShiftedLogNormalMixture(η::Array{num1},λ::Array{num2},num3) where {num1,num2,num3<: Number}
	
Where:\n
		η  =	Array of volatilities.
		λ  = 	Array of weights.
		α  = 	shift.
"""
mutable struct ShiftedLogNormalMixture{num <: Number, num2 <: Number, num3 <: Number, abstrUnderlying <: AbstractUnderlying, numtype <: Number} <: ItoProcess{numtype}
    η::Array{num, 1}
    λ::Array{num2, 1}
    α::num3
    underlying::abstrUnderlying
    function ShiftedLogNormalMixture(η::Array{num, 1}, λ::Array{num2, 1}, α::num3, underlying::abstrUnderlying) where {num <: Number, num2 <: Number, num3 <: Number, abstrUnderlying <: AbstractUnderlying}
        if minimum(η) <= 0
            error("Volatilities must be positive")
        elseif minimum(λ) <= 0
            error("weights must be positive")
        elseif sum(λ) > 1.0
            error("λs must be weights")
        elseif length(λ) != length(η) - 1
            error("Check vector lengths")
        else
            zero_typed = zero(num) + zero(num2) + zero(num3)
            return new{num, num2, num3, abstrUnderlying, typeof(zero_typed)}(η, λ, α, underlying)
        end
    end
end

export ShiftedLogNormalMixture;

function simulate!(X, mcProcess::ShiftedLogNormalMixture, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    @assert T > 0
    r = rfCurve.r
    S0 = mcProcess.underlying.S0
    d = dividend(mcProcess)
    Nstep = mcBaseData.Nstep
    η_gbm = copy(mcProcess.η)
    λ_gmb = copy(mcProcess.λ)
    mu_gbm = r - d
    dt = T / Nstep
    A0 = S0 * (1 - mcProcess.α)
    simulate!(X, LogNormalMixture(η_gbm, λ_gmb, Underlying(A0, d)), rfCurve, mcBaseData, T)
    tt = collect(0:dt:T)
    X .= X .+ mcProcess.α .* S0 .* exp.(integral(mu_gbm, t_) for t_ in tt)'
    return nothing
end
