#Parallel Modes
abstract type BaseMode end

"""
    SerialMode(seed = 0, rng = MersenneTwister())

Single threaded simulation drawing from `rng`, which is reset to `seed` by every metric.
Any `Random.AbstractRNG` can be used.
"""
struct SerialMode{rngType <: Random.AbstractRNG} <: BaseMode
    seed::Int64
    rng::rngType
    function SerialMode(seed::num = 0, rng::rngType = MersenneTwister()) where {num <: Integer, rngType <: Random.AbstractRNG}
        return new{rngType}(Int64(seed), rng)
    end
end
abstract type ParallelMode <: BaseMode end

# Defined here and not in the CUDA extension: an extension can only add methods, and the package extensions of other
# packages can dispatch on these types only if they exist without CUDA.
abstract type AbstractCudaMode <: ParallelMode end

"""
    CudaMode(seed_cpu = 0, seed_gpu = 0, rng = MersenneTwister())

Simulation on the GPU, which needs `using CUDA`. Every metric resets `rng` to `seed_cpu` and the CURAND generator to `seed_gpu`.
"""
struct CudaMode{rngType <: Random.AbstractRNG} <: AbstractCudaMode
    seed_cpu::Int64
    seed_gpu::Int64
    rng::rngType
    function CudaMode(seed_cpu::num = 0, seed_gpu::num2 = 0, rng::rngType = MersenneTwister()) where {num <: Integer, num2 <: Integer, rngType <: Random.AbstractRNG}
        return new{rngType}(Int64(seed_cpu), Int64(seed_gpu), rng)
    end
end

abstract type AbstractMultiProcess <: ParallelMode end

# Number of batches of MultiProcess by default: the workers, in the Distributed extension.
function default_nbatches end

"""
    MultiProcess(sub_mod = SerialMode(), nbatches = nworkers(), seeds = 1:nbatches)

Parallel mode splitting the paths into `nbatches` batches priced on the workers of `Distributed`, which
needs `using Distributed`, each with the parallel mode `sub_mod` and its own seed from `seeds`. The
result is the average of the batch prices.
"""
struct MultiProcess{abstractMode <: BaseMode} <: AbstractMultiProcess
    nbatches::Int64
    sub_mod::abstractMode
    seeds::Array{Int64}
    function MultiProcess(sub_mod::abstractMode = SerialMode(), nbatches::Int64 = default_nbatches(), seeds::Array{Int64} = collect(1:nbatches)) where {abstractMode <: BaseMode}
        ChainRulesCore.@ignore_derivatives @assert length(seeds) == nbatches
        return new{abstractMode}(nbatches, sub_mod, seeds)
    end
end

#Numerical Modes, keyholes for modes different than Monte Carlo
abstract type AbstractMethod end
abstract type AbstractMonteCarloMethod <: AbstractMethod end

"Plain Monte Carlo sampling (the default method of [`MonteCarloConfiguration`](@ref))."
struct StandardMC <: AbstractMonteCarloMethod end

"Antithetic variates: the second half of the paths uses the negated normals of the first half. Requires an even number of paths."
struct AntitheticMC <: AbstractMonteCarloMethod end

"Quasi Monte Carlo sampling with Sobol low discrepancy sequences."
struct SobolMode <: AbstractMonteCarloMethod end

"Reserved for simulations from prescribed random numbers; not implemented by the built-in models."
struct PrescribedMC <: AbstractMonteCarloMethod end
