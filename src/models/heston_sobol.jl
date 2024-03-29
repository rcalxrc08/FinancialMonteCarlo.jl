using Sobol, StatsFuns

function simulate!(X, mcProcess::HestonProcess, rfCurve::ZeroRate, mcBaseData::SerialSobolMonteCarloConfig, T::numb) where {numb <: Number}
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
    ChainRulesCore.@ignore_derivatives @assert T > 0.0

    ####Simulation
    ## Simulate
    κ_s = κ + λ1
    θ_s = κ * θ / κ_s
    seq = SobolSeq(2 * Nstep)
    skip(seq, Nstep * Nsim)
    dt = T / Nstep
    isDualZero = T * r * σ₀ * θ_s * κ_s * σ * ρ * 0.0 * S0
    view(X, :, 1) .= isDualZero
    vec = Array{typeof(get_rng_type(isDualZero))}(undef, 2 * Nstep)
    for i = 1:Nsim
        v_m = σ₀^2
        next!(seq, vec)
        vec .= norminvcdf.(vec)
        for j = 1:Nstep
            @views e1 = vec[j]
            @views e2 = e1 * ρ + vec[Nstep+j] * sqrt(1 - ρ * ρ)
            @views X[i, j+1] = X[i, j] + ((r - d) - 0.5 * v_m) * dt + sqrt(v_m) * sqrt(dt) * e1
            v_m += κ_s * (θ_s - v_m) * dt + σ * sqrt(v_m) * sqrt(dt) * e2
            #when v_m = 0.0, the derivative becomes NaN
            v_m = max(v_m, isDualZero)
        end
    end
    ## Conclude
    @. X = S0 * exp(X)
    return
end

function simulate!(X, mcProcess::HestonProcess, rfCurve::ZeroRateCurve, mcBaseData::SerialSobolMonteCarloConfig, T::numb) where {numb <: Number}
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
    ChainRulesCore.@ignore_derivatives @assert T > 0.0

    ####Simulation
    ## Simulate
    κ_s = κ + λ1
    θ_s = κ * θ / (κ + λ1)
    r_d = r - d
    dt = T / Nstep
    seq = SobolSeq(2 * Nstep)
    skip(seq, Nstep * Nsim)
    zero_drift = incremental_integral(r_d, dt * 0.0, dt)
    isDualZeroVol = T * σ₀ * κ * θ * λ1 * σ * ρ * 0.0 * zero_drift
    isDualZero = S0 * isDualZeroVol
    view(X, :, 1) .= isDualZero
    vec = Array{typeof(get_rng_type(isDualZero))}(undef, 2 * Nstep)
    tmp_ = [incremental_integral(r_d, (j - 1) * dt, dt) for j = 1:Nstep]
    isDualZeroVol_eps = isDualZeroVol + eps(isDualZeroVol)
    for i = 1:Nsim
        v_m = σ₀^2
        next!(seq, vec)
        @. vec = norminvcdf(vec)
        for j = 1:Nstep
            @views e1 = vec[j]
            @views e2 = e1 * ρ + vec[Nstep+j] * sqrt(1 - ρ * ρ)
            @views X[i, j+1] = X[i, j] + (tmp_[j] - 0.5 * v_m) * dt + sqrt(v_m) * sqrt(dt) * e1
            v_m += κ_s * (θ_s - v_m) * dt + σ * sqrt(v_m) * sqrt(dt) * e2
            #when v_m = 0.0, the derivative becomes NaN
            v_m = max(v_m, isDualZeroVol_eps)
        end
    end
    ## Conclude
    @. X = S0 * exp(X)
    return
end
