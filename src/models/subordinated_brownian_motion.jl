"""
Struct for SubordinatedBrownianMotion

		subordinatedBrownianMotion=SubordinatedBrownianMotion(sigma::num,drift::num1,underlying::abstrUnderlying)
	
Where:

		sigma       =	Volatility.
		drift       = 	drift.
		underlying  = 	underlying.
"""
mutable struct SubordinatedBrownianMotion{num <: Number, num1 <: Number, num2 <: Number, Distr <: Distribution{Univariate, Continuous}, numtype <: Number} <: AbstractMonteCarloProcess{numtype}
    sigma::num
    drift::num1
    θ::num2
    subordinator_::Distr
    function SubordinatedBrownianMotion(sigma::num, drift::num1, θ::num2, dist::Distr) where {num <: Number, num1 <: Number, num2 <: Number, Distr <: Distribution{Univariate, Continuous}}
        @assert sigma > 0 "volatility must be positive"
        zero_typed = zero(num) + zero(num1)
        return new{num, num1, num2, Distr, typeof(zero_typed)}(sigma, drift, θ, dist)
    end
end

export SubordinatedBrownianMotion;

function simulate!(X, mcProcess::SubordinatedBrownianMotion, mcBaseData::SerialMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    dt = T / Nstep
    drift = mcProcess.drift * dt
    sigma = mcProcess.sigma

    @assert T > 0.0

    type_sub = typeof(rand(mcBaseData.parallelMode.rng, mcProcess.subordinator_))
    isDualZero = drift * zero(type_sub)
    @views X[:, 1] .= isDualZero
    Z = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim)
    dt_s = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim)
    θ = mcProcess.θ
    @inbounds for i = 1:Nstep
        rand!(mcBaseData.parallelMode.rng, mcProcess.subordinator_, dt_s)
        randn!(mcBaseData.parallelMode.rng, Z)
        @views @. X[:, i+1] = X[:, i] + drift + θ * dt_s + sigma * sqrt(dt_s) * Z
    end
end

function simulate!(X, mcProcess::SubordinatedBrownianMotion, mcBaseData::SerialAntitheticMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    dt = T / Nstep
    drift = mcProcess.drift * dt
    sigma = mcProcess.sigma

    @assert T > 0.0

    type_sub = typeof(quantile(mcProcess.subordinator_, 0.5))
    isDualZero = drift * sigma * zero(type_sub) * 0.0
    @views X[:, 1] .= isDualZero
    Nsim_2 = div(Nsim, 2)
    Xp = @views X[1:Nsim_2, :]
    Xm = @views X[(Nsim_2+1):end, :]

    Z = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    dt_s = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    θ = mcProcess.θ
    @inbounds for j = 1:Nstep
        rand!(mcBaseData.parallelMode.rng, mcProcess.subordinator_, dt_s)
        randn!(mcBaseData.parallelMode.rng, Z)
        sqrt_dt_s = sqrt.(dt_s)
        @views @. Xp[:, j+1] = Xp[:, j] + dt_s * θ + drift + sqrt_dt_s * Z * sigma
        @views @. Xm[:, j+1] = Xm[:, j] + dt_s * θ + drift - sqrt_dt_s * Z * sigma
    end
end
