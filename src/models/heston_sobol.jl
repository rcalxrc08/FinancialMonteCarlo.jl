function simulate!(X, mcProcess::HestonProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::SobolConfig, T::Number)
    S0 = mcProcess.underlying.S0
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    ρ = mcProcess.ρ
    κ_s, θ_s, drifts, dt, isDualZeroVol = heston_setup(mcProcess, rfCurve, mcBaseData, T)
    isDualZero = isDualZeroVol * S0
    fill!(view(X, :, 1), isDualZero)
    v_floor = variance_floor(isDualZeroVol)
    ρ_orthogonal = sqrt(1 - ρ * ρ)
    sqrt_dt = sqrt(dt)
    v_0 = mcProcess.σ₀^2
    on_sobol_paths(2 * Nstep, mcBaseData, typeof(get_rng_type(isDualZero))) do i, normals
        v_m = v_0
        @inbounds for j = 1:Nstep
            e1 = normals[j]
            e2 = e1 * ρ + normals[Nstep+j] * ρ_orthogonal
            X[i, j+1] = X[i, j] + drift_of_step(drifts, j) - v_m * dt / 2 + sqrt(v_m) * sqrt_dt * e1
            v_m += κ_s * (θ_s - v_m) * dt + σ * sqrt(v_m) * sqrt_dt * e2
            v_m = max(v_m, v_floor)
        end
    end
    ## Conclude
    @. X = S0 * exp(X)
    return
end
