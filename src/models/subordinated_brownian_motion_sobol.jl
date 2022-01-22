using Sobol, StatsFuns

function simulate!(X, mcProcess::SubordinatedBrownianMotion, mcBaseData::SerialSobolMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    drift = mcProcess.drift * T / Nstep
    sigma = mcProcess.sigma
    @assert T > 0
    type_sub = typeof(quantile(mcProcess.subordinator_, 0.5))
    isDualZero = drift * zero(type_sub) * 0.0
    @views X[:, 1] .= isDualZero
    seq = SobolSeq(Nstep)
    skip(seq, Nstep * Nsim)
    vec = Array{typeof(get_rng_type(isDualZero))}(undef, Nstep)
    θ = mcProcess.θ
    @inbounds for i = 1:Nsim
        next!(seq, vec)
        @. vec = norminvcdf(vec)
        @inbounds for j = 1:Nstep
            dt_s = rand(mcBaseData.parallelMode.rng, mcProcess.subordinator_)
            @views X[i, j+1] = X[i, j] + drift + θ * dt_s + sigma * sqrt(dt_s) * vec[j]
        end
    end
end

function simulate!(X, mcProcess::SubordinatedBrownianMotionVec, mcBaseData::SerialSobolMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    drift = mcProcess.drift
    dt = T / Nstep
    sigma = mcProcess.sigma
    @assert T > 0
    type_sub = typeof(quantile(mcProcess.subordinator_, 0.5))
    zero_drift = incremental_integral(drift, zero(type_sub), zero(type_sub) + T / Nstep) * 0
    isDualZero = zero_drift * zero(type_sub) * 0
    @views X[:, 1] .= isDualZero
    seq = SobolSeq(Nstep)
    skip(seq, Nstep * Nsim)
    vec = Array{typeof(get_rng_type(isDualZero))}(undef, Nstep)
    θ = mcProcess.θ
    @inbounds for i = 1:Nsim
        # t_s = zero(type_sub)
        next!(seq, vec)
        @. vec = norminvcdf(vec)
        @inbounds for j = 1:Nstep
            dt_s = rand(mcBaseData.parallelMode.rng, mcProcess.subordinator_)
            # tmp_drift = incremental_integral(drift, t_s, dt_s)
            tmp_drift = incremental_integral(drift, (j - 1) * dt, dt)
            # t_s += dt_s
            @views X[i, j+1] = X[i, j] + tmp_drift + θ * dt_s + sigma * sqrt(dt_s) * vec[j]
        end
    end

    return X
end
