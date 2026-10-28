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

# Defined here and not with the CUDA code, that Requires loads at run time: package extensions of other packages are
# precompiled, so they can dispatch on these types only if they exist without CUDA.
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
