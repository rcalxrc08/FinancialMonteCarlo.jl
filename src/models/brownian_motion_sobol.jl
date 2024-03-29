using Sobol, StatsFuns

function simulate!(X, mcProcess::BrownianMotion, mcBaseData::SerialSobolMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    μ = mcProcess.μ
    ChainRulesCore.@ignore_derivatives @assert T > 0
    dt = T / Nstep
    mean_bm = μ * dt
    stddev_bm = σ * sqrt(dt)
    isDualZero = mean_bm * stddev_bm * 0
    view(X, :, 1) .= isDualZero
    seq = SobolSeq(Nstep)
    skip(seq, Nstep * Nsim)
    vec = Array{typeof(get_rng_type(isDualZero))}(undef, Nstep)
    @inbounds for i = 1:Nsim
        next!(seq, vec)
        @. vec = norminvcdf(vec)
        @inbounds for j = 1:Nstep
            @views X[i, j+1] = X[i, j] + mean_bm + stddev_bm * vec[j]
        end
    end
    nothing
end

function simulate!(X, mcProcess::BrownianMotionVec, mcBaseData::SerialSobolMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    μ = mcProcess.μ
    ChainRulesCore.@ignore_derivatives @assert T > 0
    dt = T / Nstep
    stddev_bm = σ * sqrt(dt)
    zero_drift = incremental_integral(μ, dt * 0, dt)
    isDualZero = stddev_bm * 0 * zero_drift
    view(X, :, 1) .= isDualZero
    seq = SobolSeq(Nstep)
    skip(seq, Nstep * Nsim)
    vec = Array{typeof(get_rng_type(isDualZero))}(undef, Nstep)
    tmp_ = [incremental_integral(μ, (j - 1) * dt, dt) for j = 1:Nstep]
    @inbounds for i = 1:Nsim
        next!(seq, vec)
        @. vec = norminvcdf(vec)
        @inbounds for j = 1:Nstep
            @views X[i, j+1] = X[i, j] + tmp_[j] + stddev_bm * vec[j]
        end
    end

    nothing
end
