
"Abstract type of the Monte Carlo pricing methods."
abstract type AbstractMonteCarloConfiguration <: AbstractMethod end

# A generator with the seed every metric resets it to.
struct SeededRNG{rngType <: Random.AbstractRNG}
    seed::Int64
    rng::rngType
end

# A generator a batch can own. The extensions of generators that cannot be copied refuse here.
batch_rng(rng::Random.AbstractRNG) = copy(rng)

# Seeded here too: with reseed = false the metrics do not, and the batches would share the same stream.
function batch_state(rng::Random.AbstractRNG, seed::Integer)
    rng_ = batch_rng(rng)
    inner_seed!(rng_, seed)
    return SeededRNG(Int64(seed), rng_)
end

# The random state of a configuration in a mode: one generator, or one per batch.
random_state(::AbstractMode, seed::Integer, rng::Random.AbstractRNG) = SeededRNG(Int64(seed), rng)
random_state(::AbstractMode, ::AbstractVector{<:Integer}, ::Random.AbstractRNG) = throw(ArgumentError("one seed per batch needs a mode with batches"))
random_state(mode::AbstractBatchedMode, seed::Integer, rng::Random.AbstractRNG) = random_state(mode, seed .+ collect(1:nbatches(mode)), rng)
function random_state(mode::AbstractBatchedMode, seeds::AbstractVector{<:Integer}, rng::Random.AbstractRNG)
    ChainRulesCore.@ignore_derivatives @assert length(seeds) == nbatches(mode) "There must be one seed per batch"
    return [batch_state(rng, seed) for seed in seeds]
end

# The seed and the generator a random state was built from (for a batched one, those of its first batch).
base_seed(state::SeededRNG) = state.seed
base_seed(states::AbstractVector{<:SeededRNG}) = first(states).seed - 1
base_rng(state::SeededRNG) = state.rng
base_rng(states::AbstractVector{<:SeededRNG}) = copy(first(states).rng)

"""
    MonteCarloConfiguration(Nsim, Nstep; sampling = StandardMC(), mode = SerialMode(), seed = 0, rng = MersenneTwister(), reseed = true)

Monte Carlo pricing method: the method argument of [`pricer`](@ref) and of the other metrics.

# Arguments
- `Nsim`: number of simulated paths, must be positive (and even for [`FinancialMonteCarlo.AntitheticMC`](@ref)).
- `Nstep`: number of time steps of each path, must be positive.
- `sampling`: an `AbstractSampling`: [`FinancialMonteCarlo.StandardMC`](@ref), [`FinancialMonteCarlo.AntitheticMC`](@ref),
  [`FinancialMonteCarlo.SobolMC`](@ref) or [`FinancialMonteCarlo.ControlVariates`](@ref).
- `mode`: an `AbstractMode`: [`FinancialMonteCarlo.SerialMode`](@ref), [`FinancialMonteCarlo.MultiThreading`](@ref),
  [`FinancialMonteCarlo.MultiProcess`](@ref) (with `Distributed`) or [`FinancialMonteCarlo.CudaMode`](@ref) (with `CUDA`).
- `seed`: the seed of `rng`. With batches (`MultiThreading`, `MultiProcess`) batch `i` has seed `seed + i`, or
  `seed[i]` when `seed` is a vector with one seed per batch.
- `rng`: the random number generator, any `Random.AbstractRNG`. Each batch draws from its own copy.
- `reseed`: when `true`, every metric resets the generators to their seeds before simulating, so the results are
  reproducible; when `false`, the streams continue from a call to the next.

# Examples
```julia
MonteCarloConfiguration(10_000, 30)
MonteCarloConfiguration(10_000, 30; seed = 42)
MonteCarloConfiguration(10_000, 30; sampling = FinancialMonteCarlo.AntitheticMC())
MonteCarloConfiguration(10_000, 30; mode = FinancialMonteCarlo.MultiThreading())
```
"""
struct MonteCarloConfiguration{samplingType <: AbstractSampling, modeType <: AbstractMode, stateType} <: AbstractMonteCarloConfiguration
    Nsim::Int64
    Nstep::Int64
    sampling::samplingType
    mode::modeType
    rng::stateType
    reseed::Bool
end

function MonteCarloConfiguration(Nsim::Integer, Nstep::Integer; sampling::AbstractSampling = StandardMC(), mode::AbstractMode = SerialMode(), seed::Union{Integer, AbstractVector{<:Integer}} = 0, rng::Random.AbstractRNG = MersenneTwister(), reseed::Bool = true)
    check_configuration(Nsim, Nstep, sampling, mode, reseed)
    return MonteCarloConfiguration(Int64(Nsim), Int64(Nstep), sampling, mode, random_state(mode, seed, rng), reseed)
end

# All the checks of a configuration, so that an unsupported combination fails when it is built.
function check_configuration(Nsim::Integer, Nstep::Integer, sampling::AbstractSampling, mode::AbstractMode, reseed::Bool)
    ChainRulesCore.@ignore_derivatives begin
        @assert Nsim > 0 "Number of Simulations must be positive"
        @assert Nstep > 0 "Number of Steps must be positive"
        check_sampling(sampling, Nsim)
        check_sampling_mode(sampling, mode)
        check_reseed(mode, reseed)
        check_mode(mode)
    end
    return nothing
end

check_sampling(::AbstractSampling, ::Integer) = nothing
check_sampling(::AntitheticMC, Nsim::Integer) = @assert iseven(Nsim) "Antithetic support only even number of simulations"

check_sampling_mode(::AbstractSampling, ::AbstractMode) = nothing
# The points of a Sobol sequence depend only on Nsim: the batches would all be the same
check_sampling_mode(::SobolMC, ::AbstractBatchedMode) = throw(ArgumentError("SobolMC does not support batches: every batch would simulate the same paths"))
check_sampling_mode(::SobolMC, ::AbstractCudaMode) = throw(ArgumentError("SobolMC is not supported in CudaMode"))

check_reseed(::AbstractMode, ::Bool) = nothing
check_reseed(::AbstractMultiProcess, reseed::Bool) = reseed || throw(ArgumentError("MultiProcess needs reseed = true: the workers receive copies of the generators, so a stream cannot continue from a call to the next"))

check_mode(::AbstractMode) = nothing
check_mode(::AbstractCudaMode) = check_extension(:FinancialMonteCarloCUDAExt, "CudaMode needs `using CUDA`")

# The generator of a configuration that is not split into batches (the batches are, see batches.jl)
host_rng(mcConfig::MonteCarloConfiguration) = mcConfig.rng.rng

function with_mode(mcConfig::MonteCarloConfiguration, mode::AbstractMode)
    return MonteCarloConfiguration(mcConfig.Nsim, mcConfig.Nstep; sampling = mcConfig.sampling, mode = mode, seed = base_seed(mcConfig.rng), rng = base_rng(mcConfig.rng), reseed = mcConfig.reseed)
end

const SobolConfig = MonteCarloConfiguration{<:SobolMC}
const BatchedConfig = MonteCarloConfiguration{<:AbstractSampling, <:AbstractBatchedMode}

export MonteCarloConfiguration;
