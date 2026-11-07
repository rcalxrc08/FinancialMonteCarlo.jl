# Names and constructors of version 0.5, to be removed in the next breaking release.
Base.@deprecate_binding SobolMode SobolMC false
Base.@deprecate_binding BaseMode AbstractMode false
Base.@deprecate_binding AbstractMonteCarloMethod AbstractSampling false

@deprecate MonteCarloConfiguration(Nsim::Integer, Nstep::Integer, seed::Integer) MonteCarloConfiguration(Nsim, Nstep; seed = seed)
@deprecate MonteCarloConfiguration(Nsim::Integer, Nstep::Integer, sampling::AbstractSampling) MonteCarloConfiguration(Nsim, Nstep; sampling = sampling)
@deprecate MonteCarloConfiguration(Nsim::Integer, Nstep::Integer, mode::AbstractMode) MonteCarloConfiguration(Nsim, Nstep; mode = mode)
@deprecate MonteCarloConfiguration(Nsim::Integer, Nstep::Integer, sampling::AbstractSampling, mode::AbstractMode) MonteCarloConfiguration(Nsim, Nstep; sampling = sampling, mode = mode)
