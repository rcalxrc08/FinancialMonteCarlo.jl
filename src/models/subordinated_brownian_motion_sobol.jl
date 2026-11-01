function simulate!(X, mcProcess::AnySubordinatedBrownianMotion, mcBaseData::SerialSobolMonteCarloConfig, T::Number)
    Nstep = mcBaseData.Nstep
    sigma = mcProcess.sigma
    θ = mcProcess.θ
    drifts, _, isDualZero = subordinated_setup(mcProcess, mcBaseData, T)
    fill!(view(X, :, 1), isDualZero)
    seq = sobol_sequence(Nstep, mcBaseData)
    vec = Array{typeof(get_rng_type(isDualZero))}(undef, Nstep)
    @inbounds for i = 1:mcBaseData.Nsim
        sobol_normals!(vec, seq)
        for j = 1:Nstep
            dt_s = rand(mcBaseData.parallelMode.rng, mcProcess.subordinator_)
            X[i, j+1] = X[i, j] + drift_of_step(drifts, j) + θ * dt_s + sigma * sqrt(dt_s) * vec[j]
        end
    end
    return
end
