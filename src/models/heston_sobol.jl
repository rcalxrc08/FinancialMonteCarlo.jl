function simulate!(X, mcProcess::HestonProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::SerialSobolMonteCarloConfig, T::Number)
    S0 = mcProcess.underlying.S0
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    ρ = mcProcess.ρ
    κ_s, θ_s, drifts, dt, isDualZeroVol = heston_setup(mcProcess, rfCurve, mcBaseData, T)
    isDualZero = isDualZeroVol * S0
    fill!(view(X, :, 1), isDualZero)
    #when v_m = 0.0, the derivative becomes NaN
    v_floor = isDualZeroVol + eps(isDualZeroVol)
    seq = sobol_sequence(2 * Nstep, mcBaseData)
    vec = Array{typeof(get_rng_type(isDualZero))}(undef, 2 * Nstep)
    ρ_orthogonal = sqrt(1 - ρ * ρ)
    sqrt_dt = sqrt(dt)
    for i = 1:mcBaseData.Nsim
        v_m = mcProcess.σ₀^2
        sobol_normals!(vec, seq)
        for j = 1:Nstep
            e1 = vec[j]
            e2 = e1 * ρ + vec[Nstep+j] * ρ_orthogonal
            X[i, j+1] = X[i, j] + drift_of_step(drifts, j) - v_m * dt / 2 + sqrt(v_m) * sqrt_dt * e1
            v_m += κ_s * (θ_s - v_m) * dt + σ * sqrt(v_m) * sqrt_dt * e2
            v_m = max(v_m, v_floor)
        end
    end
    ## Conclude
    @. X = S0 * exp(X)
    return
end
