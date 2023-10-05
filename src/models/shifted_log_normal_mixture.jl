"""
Struct for ShiftedLogNormalMixture

		lnmModel=ShiftedLogNormalMixture(η::Array{num1},λ::Array{num2},num3) where {num1,num2,num3<: Number}
	
Where:

		η  =	Array of volatilities.
		λ  = 	Array of weights.
		α  = 	shift.
"""
struct ShiftedLogNormalMixture{num <: Number, num2 <: Number, num3 <: Number, abstrUnderlying <: AbstractUnderlying, numtype <: Number} <: ItoProcess{numtype}
    η::Array{num, 1}
    λ::Array{num2, 1}
    α::num3
    underlying::abstrUnderlying
    function ShiftedLogNormalMixture(η::Array{num, 1}, λ::Array{num2, 1}, α::num3, underlying::abstrUnderlying) where {num <: Number, num2 <: Number, num3 <: Number, abstrUnderlying <: AbstractUnderlying}
        @assert minimum(η) > 0 "Volatilities must be positive"
        @assert minimum(λ) > 0 "weights must be positive"
        @assert sum(λ) <= 1.0 "λs must be weights"
        @assert length(λ) == length(η) - 1 "Check vector lengths"
        zero_typed = zero(num) + zero(num2) + zero(num3)
        return new{num, num2, num3, abstrUnderlying, typeof(zero_typed)}(η, λ, α, underlying)
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
    zero_typed = zero(typeof(integral(mu_gbm, dt) * dt))
    array_type = get_array_type(mcBaseData, mcProcess, zero_typed)
    tmp_::typeof(array_type) = exp.(integral(mu_gbm, t_) for t_ = 0:dt:T)
    @. X = X + mcProcess.α * S0 * (tmp_')
    return nothing
end
