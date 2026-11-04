function simulate!(X, mcProcess::AnySubordinatedBrownianMotion, mcBaseData::SerialSobolMonteCarloConfig, T::Number)
    Nstep = mcBaseData.Nstep
    sigma = mcProcess.sigma
    θ = mcProcess.θ
    drifts, _, isDualZero = subordinated_setup(mcProcess, mcBaseData, T)
    fill!(view(X, :, 1), isDualZero)
    rng = mcBaseData.parallelMode.rng
    subordinator_ = mcProcess.subordinator_
    on_sobol_paths(Nstep, mcBaseData, typeof(get_rng_type(isDualZero))) do i, normals
        @inbounds for j = 1:Nstep
            dt_s = rand(rng, subordinator_)
            X[i, j+1] = X[i, j] + drift_of_step(drifts, j) + θ * dt_s + sigma * sqrt(dt_s) * normals[j]
        end
    end
    return
end
