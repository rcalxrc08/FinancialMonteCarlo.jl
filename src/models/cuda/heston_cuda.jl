
function simulate!(X, mcProcess::HestonProcess, spotData::ZeroRate, mcBaseData::CudaMonteCarloConfig, T::numb) where {numb <: Number}
    r = spotData.r
    S0 = mcProcess.underlying.S0
    d = mcProcess.underlying.d
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    σ_zero = mcProcess.σ₀
    λ1 = mcProcess.λ
    κ = mcProcess.κ
    ρ = mcProcess.ρ
    θ = mcProcess.θ
    @assert T > 0

    ####Simulation
    ## Simulate
    κ_s = κ + λ1
    θ_s = κ * θ / (κ + λ1)
    r_d = r - d

    dt = T / Nstep
    isDualZeroVol = zero(T * r * σ_zero * κ * θ * λ1 * σ * ρ)
    isDualZero = zero(S0 * isDualZeroVol)
    @views fill!(X[:, 1], isDualZero) #more efficient than @views X_cu[:, 1] .= isDualZero
    v_m = CUDA.zeros(typeof(isDualZeroVol), Nsim) .+ σ_zero^2

    e1 = CuArray{typeof(get_rng_type(isDualZero))}(undef, Nsim)
    e2_rho = CuArray{typeof(get_rng_type(isDualZero))}(undef, Nsim)
    e2 = CuArray{typeof(get_rng_type(isDualZero))}(undef, Nsim)
    isDualZero_eps = isDualZeroVol + eps(isDualZeroVol)

    for j = 1:Nstep
        randn!(CUDA.CURAND.default_rng(), e1)
        randn!(CUDA.CURAND.default_rng(), e2_rho)
        @. e2 = e1 * ρ + e2_rho * sqrt(1 - ρ * ρ)
        @views @inbounds @. X[:, j+1] = X[:, j] + (r_d - 0.5 * v_m) * dt + sqrt(v_m) * sqrt(dt) * e1
        @. v_m += κ_s * (θ_s - v_m) * dt + σ * sqrt(v_m) * sqrt(dt) * e2
        #when v_m = 0.0, the derivative becomes NaN
        @. v_m = max(v_m, isDualZero_eps)
    end
    ## Conclude
    @. X = S0 * exp(X)
    return
end

function simulate!(X, mcProcess::HestonProcess, spotData::ZeroRate, mcBaseData::CudaAntitheticMonteCarloConfig, T::numb) where {numb <: Number}
    r = spotData.r
    S0 = mcProcess.underlying.S0
    d = mcProcess.underlying.d
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    σ_zero = mcProcess.σ₀
    λ1 = mcProcess.λ
    κ = mcProcess.κ
    ρ = mcProcess.ρ
    θ = mcProcess.θ
    @assert T > 0

    ####Simulation
    ## Simulate
    κ_s = κ + λ1
    θ_s = κ * θ / (κ + λ1)
    r_d = r - d

    dt = T / Nstep
    isDualZeroVol = zero(T * r * σ_zero * κ * θ * λ1 * σ * ρ)
    isDualZero = zero(S0 * isDualZeroVol)

    @views fill!(X[:, 1], isDualZero) #more efficient than @views X_cu[:, 1] .= isDualZero
    Nsim_2 = div(Nsim, 2)

    e1 = CuArray{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    e2_rho = CuArray{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    e2 = CuArray{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)

    Xp = @views X[1:Nsim_2, :]
    Xm = @views X[(Nsim_2+1):end, :]

    v_m_1 = CUDA.zeros(typeof(isDualZeroVol), Nsim_2) .+ σ_zero^2
    v_m_2 = CUDA.zeros(typeof(isDualZeroVol), Nsim_2) .+ σ_zero^2
    isDualZero_eps = isDualZeroVol + eps(isDualZeroVol)
    @inbounds for j = 1:Nstep
        randn!(CUDA.CURAND.default_rng(), e1)
        randn!(CUDA.CURAND.default_rng(), e2_rho)
        @. e2 = -(e1 * ρ + e2_rho * sqrt(1 - ρ * ρ))
        @views @. Xp[:, j+1] = Xp[:, j] + (r_d - 0.5 * v_m_1) * dt + sqrt(v_m_1) * sqrt(dt) * e1
        @views @. Xm[:, j+1] = Xm[:, j] + (r_d - 0.5 * v_m_2) * dt - sqrt(v_m_2) * sqrt(dt) * e1
        @. v_m_1 += κ_s * (θ_s - v_m_1) * dt + σ * sqrt(v_m_1) * sqrt(dt) * e2
        @. v_m_2 += κ_s * (θ_s - v_m_2) * dt - σ * sqrt(v_m_2) * sqrt(dt) * e2
        @. v_m_1 = max(v_m_1, isDualZero_eps)
        @. v_m_2 = max(v_m_2, isDualZero_eps)
    end
    ## Conclude
    @. X = S0 * exp(X)
    return
end
