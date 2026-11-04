function simulate!(X, mcProcess::Union{BrownianMotion, BrownianMotionVec}, mcBaseData::SerialSobolMonteCarloConfig, T::Number)
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    μ = mcProcess.μ
    ChainRulesCore.@ignore_derivatives @assert T > 0
    dt = T / Nstep
    stddev_bm = σ * sqrt(dt)
    drifts = step_drifts(μ, dt, Nstep)
    isDualZero = stddev_bm * 0 * drift_of_step(drifts, 1)
    fill!(view(X, :, 1), isDualZero)
    on_sobol_paths(Nstep, mcBaseData, typeof(get_rng_type(isDualZero))) do i, normals
        @inbounds for j = 1:Nstep
            X[i, j+1] = X[i, j] + drift_of_step(drifts, j) + stddev_bm * normals[j]
        end
    end
    return
end
