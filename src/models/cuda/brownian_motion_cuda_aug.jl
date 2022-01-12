
function simulate!(X_cu, mcProcess::BrownianMotionVec, mcBaseData::CudaMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep

    σ = mcProcess.σ
    μ = mcProcess.μ
    @assert T > 0
    dt = T / Nstep
    stddev_bm = σ * sqrt(dt)
    zero_drift = incremental_integral(μ, dt * 0, dt)
    isDualZero = stddev_bm * 0 * zero_drift
    Z = CuArray{typeof(get_rng_type(isDualZero))}(undef, Nsim)
    @views fill!(X_cu[:, 1], isDualZero) #more efficient than @views X_cu[:, 1] .= isDualZero
    for i = 1:Nstep
        tmp_ = incremental_integral(μ, (i - 1) * dt, dt)
        randn!(CUDA.CURAND.default_rng(), Z)
        @inbounds @. X_cu[:, i+1] = X_cu[:, i] + tmp_ + stddev_bm * Z
    end
    return
end

function simulate!(X_cu, mcProcess::BrownianMotionVec, mcBaseData::CudaAntitheticMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep

    σ = mcProcess.σ
    μ = mcProcess.μ
    @assert T > 0
    dt = T / Nstep
    stddev_bm = σ * sqrt(dt)
    zero_drift = incremental_integral(μ, dt * 0, dt)
    isDualZero = stddev_bm * 0 * zero_drift
    Nsim_2 = div(Nsim, 2)
    @views fill!(X_cu[:, 1], isDualZero) #more efficient than @views X_cu[:, 1] .= isDualZero
    X_cu_p = view(X_cu, 1:Nsim_2, :)
    X_cu_m = view(X_cu, (Nsim_2+1):Nsim, :)
    Z = CuArray{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)

    for i = 1:Nstep
        randn!(CUDA.CURAND.default_rng(), Z)
        tmp_ = incremental_integral(μ, (i - 1) * dt, dt)
        @views @. X_cu_p[:, i+1] = X_cu_p[:, i] + tmp_ + stddev_bm * Z
        @views @. X_cu_m[:, i+1] = X_cu_m[:, i] + tmp_ - stddev_bm * Z
    end
    return
end
