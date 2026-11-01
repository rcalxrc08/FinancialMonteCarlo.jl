"""
    MultiThreading(sub_mod = SerialMode(), numberOfBatch = Threads.nthreads(), seeds = sub_mod.seed .+ (1:numberOfBatch))

Parallel mode splitting the paths into `numberOfBatch` batches priced on different threads, each
with its own seed from `seeds` and its own copy of the generator of `sub_mod`, a [`FinancialMonteCarlo.SerialMode`](@ref).
The result does not depend on the number of threads. Prices are the averages of the batch prices,
weighted by their number of paths.

Generators that cannot be copied (VectorizedRNG) and [`FinancialMonteCarlo.SobolMode`](@ref) are not
supported. Of the automatic differentiation backends, only the forward ones are supported: the
batches would write to the same tape of ReverseDiff, and Mooncake cannot differentiate through tasks.
"""
struct MultiThreading{modeType <: SerialMode} <: AbstractMultiThreading
    sub_modes::Vector{modeType}
    function MultiThreading(sub_mod::SerialMode = SerialMode(), numberOfBatch::Int64 = Int64(Threads.nthreads()), seeds::Array{Int64} = sub_mod.seed .+ collect(1:numberOfBatch))
        sub_modes = batch_modes(sub_mod, numberOfBatch, seeds)
        return new{eltype(sub_modes)}(sub_modes)
    end
end

run_on_batches(f, configs, ::AbstractMultiThreading) = fetch.([Threads.@spawn f(config) for config in configs])
