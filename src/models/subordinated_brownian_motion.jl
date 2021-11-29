"""
Struct for SubordinatedBrownianMotion

		subordinatedBrownianMotion=SubordinatedBrownianMotion(sigma::num,drift::num1,underlying::abstrUnderlying)
	
Where:\n
		sigma       =	Volatility.
		drift       = 	drift.
		underlying  = 	underlying.
"""
mutable struct SubordinatedBrownianMotion{num <: Number, num1 <: Number, Distr <: Distribution{Univariate, Continuous}, numtype <: Number} <: AbstractMonteCarloProcess{numtype}
    sigma::num
    drift::num1
    subordinator_::Distr
    function SubordinatedBrownianMotion(sigma::num, drift::num1, dist::Distr) where {num <: Number, num1 <: Number, Distr <: Distribution{Univariate, Continuous}}
        @assert sigma > 0 "volatility must be positive"
        zero_typed = zero(num) + zero(num1)
        return new{num, num1, Distr, typeof(zero_typed)}(sigma, drift, dist)
    end
end

export SubordinatedBrownianMotion;

function simulate!(X, mcProcess::SubordinatedBrownianMotion, mcBaseData::MonteCarloConfiguration{type1, type2, type3, SerialMode, type5}, T::numb) where {numb <: Number, type1 <: Integer, type2 <: Integer, type3 <: StandardMC, type5 <: Random.AbstractRNG}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    drift = mcProcess.drift
    sigma = mcProcess.sigma

    @assert T > 0.0

    type_sub = typeof(rand(mcBaseData.rng, mcProcess.subordinator_))
    isDualZero = drift * zero(type_sub) * 0.0
    @views X[:, 1] .= isDualZero

    for i = 1:Nstep
        # SUBORDINATED BROWNIAN MOTION (dt_s=time change)
        for j = 1:Nsim
            dt_s = rand(mcBaseData.rng, mcProcess.subordinator_)
            @views X[j, i+1] = X[j, i] + drift * dt_s + sigma * sqrt(dt_s) * randn(mcBaseData.rng)
        end
    end

    #return X;
end

function simulate!(X, mcProcess::SubordinatedBrownianMotion, mcBaseData::MonteCarloConfiguration{type1, type2, type3, SerialMode, type5}, T::numb) where {numb <: Number, type1 <: Integer, type2 <: Integer, type3 <: AntitheticMC, type5 <: Random.AbstractRNG}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    drift = mcProcess.drift
    sigma = mcProcess.sigma

    @assert T > 0.0

    type_sub = typeof(quantile(mcProcess.subordinator_, 0.5))
    isDualZero = drift * sigma * zero(type_sub) * 0.0
    @views X[:, 1] .= isDualZero
    Nsim_2 = div(Nsim, 2)

    for j = 1:Nstep
        # SUBORDINATED BROWNIAN MOTION (dt_s=time change)
        for i = 1:Nsim_2
            dt_s = rand(mcBaseData.rng, mcProcess.subordinator_)
            Z = randn(mcBaseData.rng)
            DW = sigma * sqrt(dt_s) * Z
            muDt = drift * dt_s
            @views X[2*i-1, j+1] = X[2*i-1, j] + muDt + DW
            @views X[2*i, j+1] = X[2*i, j] + muDt - DW
        end
    end

    nothing
end
