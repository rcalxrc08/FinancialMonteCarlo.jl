"""
    MultiThreading(nbatches = Threads.nthreads())

Parallel mode splitting the paths into `nbatches` batches priced on different threads. Each batch has its own seed
and its own copy of the generator of the [`MonteCarloConfiguration`](@ref), so the result does not depend on the
number of threads. Prices are the averages of the batch prices, weighted by their number of paths.

Generators that cannot be copied (VectorizedRNG) and [`FinancialMonteCarlo.SobolMC`](@ref) are not
supported. Of the automatic differentiation backends, only the forward ones are supported: the
batches would write to the same tape of ReverseDiff, and Mooncake cannot differentiate through tasks.
"""
struct MultiThreading <: AbstractMultiThreading
    nbatches::Int64
    function MultiThreading(nbatches::Integer = Threads.nthreads())
        check_nbatches(nbatches)
        return new(Int64(nbatches))
    end
end

run_on_batches(f, configs, ::AbstractMultiThreading) = fetch.([Threads.@spawn f(config) for config in configs])
