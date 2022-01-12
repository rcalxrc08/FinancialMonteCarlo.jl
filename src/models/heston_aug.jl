
function simulate!(X, mcProcess::HestonProcess, rfCurve::ZeroRateCurve, mcBaseData::SerialMonteCarloConfig, T::numb) where {numb <: Number}
    r = rfCurve.r
    S0 = mcProcess.underlying.S0
    d = dividend(mcProcess)
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    σ₀ = mcProcess.σ₀
    λ1 = mcProcess.λ
    κ = mcProcess.κ
    ρ = mcProcess.ρ
    θ = mcProcess.θ
    @assert T > 0.0

    ## Simulate
    κ_s = κ + λ1
    θ_s = κ * θ / (κ + λ1)
    r_d = r - d
    dt = T / Nstep
    zero_drift = incremental_integral(r_d, dt * 0.0, dt)
    isDualZeroVol = zero(T * σ₀ * κ * θ * λ1 * σ * ρ * zero_drift)
    isDualZero = S0 * isDualZeroVol
    X[:, 1] .= isDualZero
    v_m = [σ₀^2 + isDualZero for _ = 1:Nsim]
    isDualZeroVol_eps = isDualZeroVol + eps(isDualZeroVol)
    e1 = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim)
    e2_rho = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim)
    e2 = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim)
    @inbounds for j = 1:Nstep
        randn!(mcBaseData.parallelMode.rng, e1)
        randn!(mcBaseData.parallelMode.rng, e2_rho)
        @. e2 = e1 * ρ + e2_rho * sqrt(1 - ρ * ρ)
        tmp_ = incremental_integral(r_d, (j - 1) * dt, dt)
        @. @views X[:, j+1] = X[:, j] + (tmp_ - 0.5 * v_m * dt) + sqrt(v_m) * sqrt(dt) * e1
        @. v_m += κ_s * (θ_s - v_m) * dt + σ * sqrt(v_m) * sqrt(dt) * e2
        @. v_m = max(v_m, isDualZeroVol_eps)
    end
    ## Conclude
    @. X = S0 * exp(X)
    return
end

function simulate!(X, mcProcess::HestonProcess, rfCurve::ZeroRateCurve, mcBaseData::SerialAntitheticMonteCarloConfig, T::numb) where {numb <: Number}
    r = rfCurve.r
    S0 = mcProcess.underlying.S0
    d = dividend(mcProcess)
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    σ₀ = mcProcess.σ₀
    λ1 = mcProcess.λ
    κ = mcProcess.κ
    ρ = mcProcess.ρ
    θ = mcProcess.θ
    @assert T > 0.0

    ####Simulation
    ## Simulate
    κ_s = κ + λ1
    θ_s = κ * θ / (κ + λ1)

    dt = T / Nstep
    isDualZeroVol = zero(T * σ₀ * κ * θ * λ1 * σ * ρ)
    isDualZero = S0 * isDualZeroVol

    X[:, 1] .= isDualZero
    Nsim_2 = div(Nsim, 2)
    r_d = r - d
    Xp = @views X[1:Nsim_2, :]
    Xm = @views X[(Nsim_2+1):end, :]
    v_mp = [σ₀^2 + isDualZero for _ = 1:Nsim_2]
    v_mm = [σ₀^2 + isDualZero for _ = 1:Nsim_2]
    e1 = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    e2_rho = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    e2 = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    isDualZeroVol_eps = isDualZeroVol + eps(isDualZeroVol)
    for j = 1:Nstep
        randn!(mcBaseData.parallelMode.rng, e1)
        randn!(mcBaseData.parallelMode.rng, e2_rho)
        @. e2 = e1 * ρ + e2_rho * sqrt(1 - ρ * ρ)
        tmp_ = incremental_integral(r_d, (j - 1) * dt, dt)
        @. @views Xp[:, j+1] = Xp[:, j] + (tmp_ - 0.5 * v_mp * dt) + sqrt(v_mp) * sqrt(dt) * e1
        @. @views Xm[:, j+1] = Xm[:, j] + (tmp_ - 0.5 * v_mm * dt) - sqrt(v_mm) * sqrt(dt) * e1
        @. v_mp += κ_s * (θ_s - v_mp) * dt + (σ * sqrt(dt)) * sqrt(v_mp) * e2
        @. v_mm += κ_s * (θ_s - v_mm) * dt - (σ * sqrt(dt)) * sqrt(v_mm) * e2
        @. v_mp = max(v_mp, isDualZeroVol_eps)
        @. v_mm = max(v_mm, isDualZeroVol_eps)
    end

    ## Conclude
    @. X = S0 * exp(X)
    return
end
