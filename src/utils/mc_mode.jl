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

# The modes splitting the paths into batches, see batches.jl
abstract type AbstractBatchedMode <: ParallelMode end
abstract type AbstractMultiThreading <: AbstractBatchedMode end
abstract type AbstractMultiProcess <: AbstractBatchedMode end

# Number of batches of MultiProcess by default: the workers, in the Distributed extension.
function default_nbatches end

function check_distributed()
    ChainRulesCore.@ignore_derivatives isnothing(Base.get_extension(@__MODULE__, :FinancialMonteCarloDistributedExt)) && throw(ArgumentError("MultiProcess needs `using Distributed`"))
end

"""
    MultiProcess(sub_mod = SerialMode(), nbatches = nworkers(), seeds = sub_mod.seed .+ (1:nbatches))

Parallel mode splitting the paths into `nbatches` batches priced on the workers of `Distributed`, which
needs `using Distributed`, each with its own seed from `seeds` and its own copy of the generator of
`sub_mod`, a [`FinancialMonteCarlo.SerialMode`](@ref). The result does not depend on the number of workers, which take the
next batch when they finish one. Prices are the averages of the batch prices, weighted by their number of paths.

The workers must load the package and everything the priced objects use (`@everywhere using FinancialMonteCarlo`).
To differentiate through a pricer, the function must exist on the workers too (`@everywhere f(x) = ...`):
an anonymous function defined only on the main process raises an error.
`simulate` and `distribution` send all the paths back to the main process, which can take longer
than simulating them; the other metrics send back only their results.

`init_rng = false`, generators that cannot be copied (VectorizedRNG) and [`FinancialMonteCarlo.SobolMode`](@ref)
are not supported. Of the automatic differentiation backends, only the forward ones are supported.
"""
struct MultiProcess{modeType <: SerialMode} <: AbstractMultiProcess
    sub_modes::Vector{modeType}
    function MultiProcess(sub_mod::SerialMode = SerialMode(), nbatches::Union{Int64, Nothing} = nothing, seeds::Union{Array{Int64}, Nothing} = nothing)
        check_distributed()
        nbatches_ = something(nbatches, default_nbatches())
        sub_modes = batch_modes(sub_mod, nbatches_, something(seeds, sub_mod.seed .+ collect(1:nbatches_)))
        return new{eltype(sub_modes)}(sub_modes)
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
