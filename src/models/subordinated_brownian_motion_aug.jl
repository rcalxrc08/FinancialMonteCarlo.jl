"""
Struct for SubordinatedBrownianMotion

		subordinatedBrownianMotion=SubordinatedBrownianMotion(sigma::num,drift::num1)
	
Where:

		sigma       =	Volatility.
		drift       = 	drift.
		underlying  = 	underlying.
"""
mutable struct SubordinatedBrownianMotionVec{num <: Number, num1 <: Number, num4 <: Number, num2 <: Number, Distr <: Distribution{Univariate, Continuous}, numtype <: Number} <: AbstractMonteCarloProcess{numtype}
    sigma::num
    drift::CurveType{num1, num4}
    θ::num2
    subordinator_::Distr
    function SubordinatedBrownianMotionVec(sigma::num, drift::CurveType{num1, num4}, θ::num2, dist::Distr) where {num <: Number, num1 <: Number, num2 <: Number, num4 <: Number, Distr <: Distribution{Univariate, Continuous}}
        @assert sigma > 0 "volatility must be positive"
        zero_typed = zero(num) + zero(num1) + zero(num4)
        return new{num, num1, num4, num2, Distr, typeof(zero_typed)}(sigma, drift, θ, dist)
    end
end

function SubordinatedBrownianMotion(σ::num, drift::CurveType{num1, num4}, θ::num2, subordinator_::Distr) where {num <: Number, num1 <: Number, num2 <: Number, num4 <: Number, Distr <: Distribution{Univariate, Continuous}}
    return SubordinatedBrownianMotionVec(σ, drift, θ, subordinator_)
end

function simulate!(X, mcProcess::SubordinatedBrownianMotionVec, mcBaseData::SerialMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep

    drift = mcProcess.drift
    sigma = mcProcess.sigma
    @assert T > 0
    type_sub = typeof(quantile(mcProcess.subordinator_, 0.5))
    dt_s = Array{type_sub}(undef, Nsim)
    dt = T / Nstep
    zero_drift = incremental_integral(drift, zero(type_sub), zero(type_sub) + dt) * 0
    isDualZero = sigma * zero(type_sub) * zero_drift
    @views X[:, 1] .= isDualZero
    Z = Array{Float64}(undef, Nsim)
    θ = mcProcess.θ
    @inbounds for i = 1:Nstep
        randn!(mcBaseData.parallelMode.rng, Z)
        rand!(mcBaseData.parallelMode.rng, mcProcess.subordinator_, dt_s)
        tmp_drift = incremental_integral(drift, (i - 1) * dt, dt)
        # SUBORDINATED BROWNIAN MOTION (dt_s=time change)
        @views @. X[:, i+1] = X[:, i] + tmp_drift + θ * dt_s + sigma * sqrt(dt_s) * Z
    end

    return
end

function simulate!(X, mcProcess::SubordinatedBrownianMotionVec, mcBaseData::SerialAntitheticMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    type_sub = typeof(quantile(mcProcess.subordinator_, 0.5))
    drift = mcProcess.drift
    sigma = mcProcess.sigma
    @assert T > 0
    dt = T / Nstep
    zero_drift = incremental_integral(drift, zero(type_sub), zero(type_sub) + dt) * 0
    isDualZero = sigma * zero(type_sub) * zero_drift
    @views X[:, 1] .= isDualZero
    Nsim_2 = div(Nsim, 2)
    dt_s = Array{type_sub}(undef, Nsim_2)
    Xp = @views X[1:Nsim_2, :]
    Xm = @views X[(Nsim_2+1):end, :]
    Z = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    tmp_drift = Array{typeof(zero_drift)}(undef, Nsim_2)
    θ = mcProcess.θ
    @inbounds for i = 1:Nstep
        randn!(mcBaseData.parallelMode.rng, Z)
        rand!(mcBaseData.parallelMode.rng, mcProcess.subordinator_, dt_s)
        # tmp_drift .= [incremental_integral(drift, t_, dt_) for (t_, dt_) in zip(t_s, dt_s)]
        tmp_drift = incremental_integral(drift, (i - 1) * dt, dt)
        # @. t_s += dt_s
        # SUBORDINATED BROWNIAN MOTION (dt_s=time change)
        @. @views Xp[:, i+1] = Xp[:, i] + tmp_drift + dt_s * θ + sigma * sqrt(dt_s) * Z
        @. @views Xm[:, i+1] = Xm[:, i] + tmp_drift + dt_s * θ - sigma * sqrt(dt_s) * Z
    end

    return
end
