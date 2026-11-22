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
            # the scalar steps of the Monte Carlo kernel (heston.jl): one formula and one rounding with it and its adjoint
            X[i, j+1] = heston_log_step(X[i, j], v_m, e1, drift_of_step(drifts, j), dt, sqrt_dt)
            v_m = max(heston_variance_step(v_m, e1, normals[Nstep+j], dt, sqrt_dt, κ_s, θ_s, σ, ρ, ρ_orthogonal), v_floor)
        end
    end
    ## Conclude
    exp_paths!(X, S0)
    return
end
