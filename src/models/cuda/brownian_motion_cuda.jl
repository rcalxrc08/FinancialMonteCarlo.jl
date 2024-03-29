
function simulate!(X_cu, mcProcess::BrownianMotion, mcBaseData::CudaMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    μ = mcProcess.μ
    ChainRulesCore.@ignore_derivatives @assert T > 0
    dt = T / Nstep
    mean_bm = μ * dt
    stddev_bm = σ * sqrt(dt)
    isDualZero = mean_bm * stddev_bm * zero(Int8)
    @views fill!(X_cu[:, 1], isDualZero) #more efficient than @views X_cu[:, 1] .= isDualZero
    Z = CuArray{typeof(get_rng_type(isDualZero))}(undef, Nsim)
    @inbounds for i = 1:Nstep
        randn!(CUDA.CURAND.default_rng(), Z)
        @views @. X_cu[:, i+1] = X_cu[:, i] + mean_bm + stddev_bm * Z
    end
    return
end

function simulate!(X_cu, mcProcess::BrownianMotion, mcBaseData::CudaAntitheticMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep

    σ = mcProcess.σ
    μ = mcProcess.μ
    ChainRulesCore.@ignore_derivatives @assert T > 0
    dt = T / Nstep
    mean_bm = μ * dt
    stddev_bm = σ * sqrt(dt)
    isDualZero = mean_bm * stddev_bm * zero(Int8)
    Nsim_2 = div(Nsim, 2)
    @views fill!(X_cu[:, 1], isDualZero) #more efficient than @views X_cu[:, 1] .= isDualZero
    X_cu_p = view(X_cu, 1:Nsim_2, :)
    X_cu_m = view(X_cu, (Nsim_2+1):Nsim, :)
    Z = CuArray{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    for i = 1:Nstep
        randn!(CUDA.CURAND.default_rng(), Z)
        @views @inbounds @. X_cu_p[:, i+1] = X_cu_p[:, i] + mean_bm + stddev_bm * Z
        @views @inbounds @. X_cu_m[:, i+1] = X_cu_m[:, i] + mean_bm - stddev_bm * Z
    end
    return
end
