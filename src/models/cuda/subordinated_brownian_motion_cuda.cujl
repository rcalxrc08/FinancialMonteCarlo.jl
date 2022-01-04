
function simulate!(X_cu, mcProcess::SubordinatedBrownianMotion, mcBaseData::MonteCarloConfiguration{type1, type2, type3, CudaMode, type5}, T::numb) where {numb <: Number, type1 <: Number, type2 <: Number, type3 <: StandardMC, type5 <: Random.AbstractRNG}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    drift = mcProcess.drift
    sigma = mcProcess.sigma
    @assert T > 0
    type_sub = typeof(Float32(quantile(mcProcess.subordinator_, 0.5)))
    dt_s = CuArray{type_sub}(undef, Nsim)
    dt_s_cpu = Array{type_sub}(undef, Nsim)
    isDualZero = drift * sigma * zero(eltype(dt_s)) * 0
    @views fill!(X_cu[:, 1], isDualZero) #more efficient than @views X_cu[:, 1] .= isDualZero
    Z = CuArray{typeof(get_rng_type(isDualZero))}(undef, Nsim)
    @inbounds for i = 1:Nstep
        # SUBORDINATED BROWNIAN MOTION (dt_s=time change)
        rand!(mcBaseData.rng, mcProcess.subordinator_, dt_s_cpu)
        dt_s .= cu(dt_s_cpu)
        randn!(CUDA.CURAND.default_rng(), Z)
        @views @. X_cu[:, i+1] = X_cu[:, i] + drift * dt_s + sigma * sqrt(dt_s) * Z
    end
    return
end

function simulate!(X_cu, mcProcess::SubordinatedBrownianMotion, mcBaseData::MonteCarloConfiguration{type1, type2, type3, CudaMode, type5}, T::numb) where {numb <: Number, type1 <: Number, type2 <: Number, type3 <: AntitheticMC, type5 <: Random.AbstractRNG}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    drift = mcProcess.drift
    sigma = mcProcess.sigma
    @assert T > 0
    type_sub = typeof(Float32(quantile(mcProcess.subordinator_, 0.5)))
    Nsim_2 = div(Nsim, 2)
    dt_s = CuArray{type_sub}(undef, Nsim_2)
    dt_s_cpu = Array{type_sub}(undef, Nsim_2)
    isDualZero = drift * sigma * zero(eltype(dt_s)) * 0
    @views fill!(X_cu[:, 1], isDualZero) #more efficient than @views X_cu[:, 1] .= isDualZero
    X_cu_p = view(X_cu, 1:Nsim_2, :)
    X_cu_m = view(X_cu, (Nsim_2+1):Nsim, :)
    Z = CuArray{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    @inbounds for i = 1:Nstep
        randn!(CUDA.CURAND.default_rng(), Z)
        # SUBORDINATED BROWNIAN MOTION (dt_s=time change)
        rand!(mcBaseData.rng, mcProcess.subordinator_, dt_s_cpu)
        dt_s .= cu(dt_s_cpu)
        @views @. X_cu_p[:, i+1] = X_cu_p[:, i] + drift * dt_s + sigma * sqrt(dt_s) * Z
        @views @. X_cu_m[:, i+1] = X_cu_m[:, i] + drift * dt_s - sigma * sqrt(dt_s) * Z
    end
    return
end
