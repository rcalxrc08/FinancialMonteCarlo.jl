"""
Struct for SubordinatedBrownianMotion

		subordinatedBrownianMotion=SubordinatedBrownianMotion(sigma::num,drift::num1)
	
Where:\n
		sigma       =	Volatility.
		drift       = 	drift.
		underlying  = 	underlying.
"""
mutable struct SubordinatedBrownianMotionVec{num <: Number, num1 <: Number, num4 <: Number, Distr <: Distribution{Univariate, Continuous}, numtype <: Number} <: AbstractMonteCarloProcess{numtype}
    sigma::num
    drift::CurveType{num1, num4}
    subordinator_::Distr
    function SubordinatedBrownianMotionVec(sigma::num, drift::CurveType{num1, num4}, dist::Distr) where {num <: Number, num1 <: Number, num4 <: Number, Distr <: Distribution{Univariate, Continuous}}
        @assert sigma > 0 "volatility must be positive"
        zero_typed = zero(num) + zero(num1) + zero(num4)
        return new{num, num1, num4, Distr, typeof(zero_typed)}(sigma, drift, dist)
    end
end

function SubordinatedBrownianMotion(σ::num, drift::CurveType{num1, num4}, subordinator_::Distr) where {num <: Number, num1 <: Number, num4 <: Number, Distr <: Distribution{Univariate, Continuous}}
    return SubordinatedBrownianMotionVec(σ, drift, subordinator_)
end

function simulate!(X, mcProcess::SubordinatedBrownianMotionVec, mcBaseData::MonteCarloConfiguration{type1, type2, type3, SerialMode, type5}, T::numb) where {numb <: Number, type1 <: Number, type2 <: Number, type3 <: StandardMC, type5 <: Random.AbstractRNG}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep

    drift = mcProcess.drift
    sigma = mcProcess.sigma
    @assert T > 0
    type_sub = typeof(quantile(mcProcess.subordinator_, 0.5))
    dt_s = Array{type_sub}(undef, Nsim)
    t_s = zeros(type_sub, Nsim)
    dt = T / Nstep
    zero_drift = incremental_integral(drift, zero(type_sub), zero(type_sub) + dt) * 0
    isDualZero = sigma * zero(type_sub) * 0 * zero_drift
    @views X[:, 1] .= isDualZero
    Z = Array{Float64}(undef, Nsim)
    @inbounds for i = 1:Nstep
        randn!(mcBaseData.rng, Z)
        rand!(mcBaseData.rng, mcProcess.subordinator_, dt_s)
        tmp_drift = [incremental_integral(drift, t_, dt_) for (t_, dt_) in zip(t_s, dt_s)]
        @. t_s += dt_s
        # SUBORDINATED BROWNIAN MOTION (dt_s=time change)
        @views @. X[:, i+1] = X[:, i] + tmp_drift + sigma * sqrt(dt_s) * Z
    end

    return
end

function simulate!(X, mcProcess::SubordinatedBrownianMotionVec, mcBaseData::MonteCarloConfiguration{type1, type2, type3, SerialMode, type5}, T::numb) where {numb <: Number, type1 <: Number, type2 <: Number, type3 <: AntitheticMC, type5 <: Random.AbstractRNG}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    type_sub = typeof(quantile(mcProcess.subordinator_, 0.5))
    drift = mcProcess.drift
    sigma = mcProcess.sigma
    @assert T > 0
    dt = T / Nstep
    zero_drift = incremental_integral(drift, zero(type_sub), zero(type_sub) + dt)
    isDualZero = sigma * zero(type_sub) * 0 * zero_drift
    #X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
    @views X[:, 1] .= isDualZero
    Nsim_2 = div(Nsim, 2)
    dt_s = Array{type_sub}(undef, Nsim_2)
    t_s = zeros(type_sub, Nsim_2)
    Xp = @views X[1:Nsim_2, :]
    Xm = @views X[(Nsim_2+1):end, :]
    Z = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    @inbounds for i = 1:Nstep
        randn!(mcBaseData.rng, Z)
        rand!(mcBaseData.rng, mcProcess.subordinator_, dt_s)
        tmp_drift = [incremental_integral(drift, t_, dt_) for (t_, dt_) in zip(t_s, dt_s)]
        @. t_s += dt_s
        # SUBORDINATED BROWNIAN MOTION (dt_s=time change)
        @. @views Xp[:, i+1] = Xp[:, i] + tmp_drift + sigma * sqrt(dt_s) * Z
        @. @views Xm[:, i+1] = Xm[:, i] + tmp_drift - sigma * sqrt(dt_s) * Z
    end

    return
end
