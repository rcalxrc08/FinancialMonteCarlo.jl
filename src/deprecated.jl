# Names and constructors of version 0.5, to be removed in the next breaking release.
Base.@deprecate_binding SobolMode SobolMC false
Base.@deprecate_binding BaseMode AbstractMode false
Base.@deprecate_binding AbstractMonteCarloMethod AbstractSampling false

@deprecate MonteCarloConfiguration(Nsim::Integer, Nstep::Integer, seed::Integer) MonteCarloConfiguration(Nsim, Nstep; seed = seed)
@deprecate MonteCarloConfiguration(Nsim::Integer, Nstep::Integer, sampling::AbstractSampling) MonteCarloConfiguration(Nsim, Nstep; sampling = sampling)
@deprecate MonteCarloConfiguration(Nsim::Integer, Nstep::Integer, mode::AbstractMode) MonteCarloConfiguration(Nsim, Nstep; mode = mode)
@deprecate MonteCarloConfiguration(Nsim::Integer, Nstep::Integer, sampling::AbstractSampling, mode::AbstractMode) MonteCarloConfiguration(Nsim, Nstep; sampling = sampling, mode = mode)

# Deprecated in 0.7: the shifted lognormal mixture is ShiftedProcess around a LogNormalMixture. The old name is an alias
# of those types, so `x isa ShiftedLogNormalMixture` still works, and its constructor builds one.
function (::Type{ShiftedProcess{<:LogNormalMixture}})(η::AbstractVector, λ::AbstractVector, α::Number, underlying::AbstractUnderlying)
    # the deprecated binding alone warns only where the name is resolved, and not at all on recent Julia versions
    Base.depwarn("`ShiftedLogNormalMixture(η, λ, α, underlying)` is deprecated, use `ShiftedProcess(LogNormalMixture(η, λ, underlying), α)` instead.", :ShiftedLogNormalMixture)
    return ShiftedProcess(LogNormalMixture(η, λ, underlying), α)
end
Base.@deprecate_binding ShiftedLogNormalMixture ShiftedProcess{<:LogNormalMixture}
