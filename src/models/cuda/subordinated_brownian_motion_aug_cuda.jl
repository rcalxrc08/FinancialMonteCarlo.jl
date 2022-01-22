
function simulate!(X, mcProcess::SubordinatedBrownianMotionVec, mcBaseData::CudaMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep

    drift = mcProcess.drift
    sigma = mcProcess.sigma
    @assert T > 0
    type_sub = typeof(quantile(mcProcess.subordinator_, 0.5))
    dt_s_cpu = Array{type_sub}(undef, Nsim)
    dt_s = CuArray{type_sub}(undef, Nsim)
    # t_s = zeros(type_sub, Nsim)
    dt = T / Nstep
    zero_drift = incremental_integral(drift, zero(type_sub), zero(type_sub) + dt) * 0
    isDualZero = sigma * zero(type_sub) * 0 * zero_drift
    @views fill!(X[:, 1], isDualZero) #more efficient than @views X[:, 1] .= isDualZero
    Z = CuArray{typeof(get_rng_type(isDualZero))}(undef, Nsim)
    # tmp_drift = Array{typeof(zero_drift)}(undef, Nsim)
    # tmp_drift_gpu = CuArray{typeof(zero_drift)}(undef, Nsim)
    θ = mcProcess.θ
    @inbounds for i = 1:Nstep
        randn!(CUDA.CURAND.default_rng(), Z)
        rand!(mcBaseData.parallelMode.rng, mcProcess.subordinator_, dt_s_cpu)
        # tmp_drift .= [incremental_integral(drift, t_, dt_) for (t_, dt_) in zip(t_s, dt_s_cpu)]
        tmp_drift = incremental_integral(drift, (i - 1) * dt, dt)
        # tmp_drift_gpu .= cu(tmp_drift)
        # @. t_s += dt_s_cpu
        dt_s .= cu(dt_s_cpu)
        # SUBORDINATED BROWNIAN MOTION (dt_s=time change)
        @views @. X[:, i+1] = X[:, i] + tmp_drift + θ * dt_s + sigma * sqrt(dt_s) * Z
    end

    return
end

function simulate!(X, mcProcess::SubordinatedBrownianMotionVec, mcBaseData::CudaAntitheticMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    type_sub = typeof(quantile(mcProcess.subordinator_, 0.5))
    drift = mcProcess.drift
    sigma = mcProcess.sigma
    @assert T > 0
    dt = T / Nstep
    zero_drift = incremental_integral(drift, zero(type_sub), zero(type_sub) + dt)
    isDualZero = sigma * zero(type_sub) * 0 * zero_drift
    @views fill!(X[:, 1], isDualZero) #more efficient than @views X[:, 1] .= isDualZero
    Nsim_2 = div(Nsim, 2)
    dt_s = CuArray{type_sub}(undef, Nsim_2)
    dt_s_cpu = Array{type_sub}(undef, Nsim_2)
    # tmp_drift = Array{typeof(zero_drift)}(undef, Nsim_2)
    # tmp_drift_gpu = CuArray{typeof(zero_drift)}(undef, Nsim_2)
    # t_s = zeros(type_sub, Nsim_2)
    Xp = @views X[1:Nsim_2, :]
    Xm = @views X[(Nsim_2+1):end, :]
    Z = CuArray{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    θ = mcProcess.θ
    @inbounds for i = 1:Nstep
        randn!(CUDA.CURAND.default_rng(), Z)
        rand!(mcBaseData.parallelMode.rng, mcProcess.subordinator_, dt_s_cpu)
        # tmp_drift .= [incremental_integral(drift, t_, dt_) for (t_, dt_) in zip(t_s, dt_s_cpu)]
        tmp_drift = incremental_integral(drift, (i - 1) * dt, dt)
        # @. t_s += dt_s_cpu
        dt_s .= cu(dt_s_cpu)
        # tmp_drift_gpu .= cu(tmp_drift)
        # SUBORDINATED BROWNIAN MOTION (dt_s=time change)
        @. @views Xp[:, i+1] = Xp[:, i] + tmp_drift + θ * dt_s + sigma * sqrt(dt_s) * Z
        @. @views Xm[:, i+1] = Xm[:, i] + tmp_drift + θ * dt_s - sigma * sqrt(dt_s) * Z
    end

    return
end
