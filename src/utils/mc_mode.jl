# Pricing methods, execution modes and sampling schemes.
#
# - A method (AbstractMethod) is what pricer dispatches on: a MonteCarloConfiguration here, a Fourier method in
#   FinancialFFT, a Vibrato method in VibratoMonteCarlo. Every method carries its execution mode.
# - A mode (AbstractMode) says where and how a method runs: serial, on threads or workers, on the GPU. It holds no
#   seed and no generator, so the methods without randomness can use it too.
# - A sampling scheme (AbstractSampling) says how a Monte Carlo method draws its paths.

"Abstract type of the pricing methods: every package implementing `pricer(model, zero_rate, method, payoff)` subtypes it."
abstract type AbstractMethod end

"""
    execution_mode(method)

The [`FinancialMonteCarlo.AbstractMode`](@ref) the pricing `method` runs in.
"""
execution_mode(method::AbstractMethod) = method.mode

"""
    with_mode(method, mode)

A copy of the pricing `method` that runs in `mode`.
"""
with_mode(method::AbstractMethod, mode) = @set method.mode = mode

## Execution modes
"Abstract type of the execution modes: where and how a pricing method runs."
abstract type AbstractMode end

"""
    SerialMode()

Execution on one thread of the CPU (the default mode of every method).
"""
struct SerialMode <: AbstractMode end

abstract type ParallelMode <: AbstractMode end

# Defined here and not in the CUDA extension: an extension can only add methods, and the package extensions of other
# packages can dispatch on these types only if they exist without CUDA.
abstract type AbstractCudaMode <: ParallelMode end

"""
    CudaMode()

Execution on the GPU, which needs `using CUDA`. The paths, and the results of a smile, stay on the device.
"""
struct CudaMode <: AbstractCudaMode end

# Where the arrays of a mode live. An input is moved there when a metric starts, once (to_mode): a smile with its strikes
# on the host is priced in CudaMode, and one with its strikes on the GPU in the other modes. The results stay where the
# mode computes them. Adapt.adapt rebuilds the structs whose arrays move; the CUDA extension adds the adaptors.
struct HostAdaptor end
mode_adaptor(::AbstractMode) = HostAdaptor()
to_mode(x, mode::AbstractMode) = Adapt.adapt(mode_adaptor(mode), x)

# The modes splitting the paths into batches, see batches.jl
abstract type AbstractBatchedMode <: ParallelMode end
abstract type AbstractMultiThreading <: AbstractBatchedMode end
abstract type AbstractMultiProcess <: AbstractBatchedMode end

nbatches(mode::AbstractBatchedMode) = mode.nbatches

function check_nbatches(nbatches::Integer)
    ChainRulesCore.@ignore_derivatives @assert nbatches > 0 "Number of batches must be positive"
end

# Number of batches of MultiProcess by default: the workers, in the Distributed extension.
function default_nbatches end

check_extension(name::Symbol, message) = ChainRulesCore.@ignore_derivatives isnothing(Base.get_extension(@__MODULE__, name)) && throw(ArgumentError(message))
check_distributed() = check_extension(:FinancialMonteCarloDistributedExt, "MultiProcess needs `using Distributed`")

"""
    MultiProcess(nbatches = nworkers())

Parallel mode splitting the paths into `nbatches` batches priced on the workers of `Distributed`, which
needs `using Distributed`. Each batch has its own seed and its own copy of the generator of the
[`MonteCarloConfiguration`](@ref), so the result does not depend on the number of workers, which take the
next batch when they finish one. Prices are the averages of the batch prices, weighted by their number of paths.

The workers must load the package and everything the priced objects use (`@everywhere using FinancialMonteCarlo`).
To differentiate through a pricer, the function must exist on the workers too (`@everywhere f(x) = ...`):
an anonymous function defined only on the main process raises an error.
`simulate` and `distribution` send all the paths back to the main process, which can take longer
than simulating them; the other metrics send back only their results.

`reseed = false`, generators that cannot be copied (VectorizedRNG) and [`FinancialMonteCarlo.SobolMC`](@ref)
are not supported. Of the automatic differentiation backends, only the forward ones are supported.
"""
struct MultiProcess <: AbstractMultiProcess
    nbatches::Int64
    function MultiProcess(nbatches::Union{Integer, Nothing} = nothing)
        check_distributed()
        nbatches_ = Int64(something(nbatches, default_nbatches()))
        check_nbatches(nbatches_)
        return new(nbatches_)
    end
end

## Sampling schemes of the Monte Carlo methods
"Abstract type of the sampling schemes of a [`MonteCarloConfiguration`](@ref)."
abstract type AbstractSampling end

"Plain Monte Carlo sampling (the default sampling of [`MonteCarloConfiguration`](@ref))."
struct StandardMC <: AbstractSampling end

"Antithetic variates: the second half of the paths uses the negated normals of the first half. Requires an even number of paths."
struct AntitheticMC <: AbstractSampling end

"Quasi Monte Carlo sampling with Sobol low discrepancy sequences. Only in [`FinancialMonteCarlo.SerialMode`](@ref)."
struct SobolMC <: AbstractSampling end
