abstract type AbstractMultiThreading <: ParallelMode end
struct MultiThreading{abstractMode <: BaseMode} <: AbstractMultiThreading
    numberOfBatch::Int64
    sub_mod::abstractMode
    seeds::Array{Int64}
    function MultiThreading(sub_mod::abstractMode = SerialMode(), numberOfBatch::Int64 = Int64(Threads.nthreads()), seeds::Array{Int64} = collect(1:Threads.nthreads())) where {abstractMode <: BaseMode}
        ChainRulesCore.@ignore_derivatives @assert length(seeds) == numberOfBatch
        return new{abstractMode}(numberOfBatch, sub_mod, seeds)
    end
end

function pricer_macro_multithreading(model_type, payoff_type)
    @eval begin
        function pricer(mcProcess::$model_type, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:AbstractMultiThreading}, abstractPayoff::$payoff_type)
            zero_typed = predict_output_type_zero(mcProcess, rfCurve, mcConfig, abstractPayoff)
            numberOfBatch = mcConfig.parallelMode.numberOfBatch
            price = zeros(typeof(zero_typed), numberOfBatch)
            Threads.@threads for i = 1:numberOfBatch
                mc_config_i = MonteCarloConfiguration(div(mcConfig.Nsim, numberOfBatch), mcConfig.Nstep, mcConfig.monteCarloMethod, mcConfig.parallelMode.sub_mod)
                mc_config_i = @views @set mc_config_i.parallelMode.seed = mcConfig.parallelMode.seeds[i]
                @views price[i] = pricer(mcProcess, rfCurve, mc_config_i, abstractPayoff)
            end
            Out = sum(price) / numberOfBatch
            return Out
        end
    end
end

pricer_macro_multithreading(BaseProcess, AbstractPayoff)
pricer_macro_multithreading(BaseProcess, Dict{AbstractPayoff, Number})
pricer_macro_multithreading(Dict{String, AbstractMonteCarloProcess}, Dict{String, Dict{AbstractPayoff, Number}})
pricer_macro_multithreading(VectorialMonteCarloProcess, Array{Dict{AbstractPayoff, Number}})

function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:AbstractMultiThreading}, abstractPayoffs::Array{abstractPayoff_}) where {abstractPayoff_ <: AbstractPayoff}
    zero_typed = predict_output_type_zero(mcProcess, rfCurve, mcConfig, abstractPayoffs)
    numberOfBatch = mcConfig.parallelMode.numberOfBatch
    price = zeros(typeof(zero_typed), length(abstractPayoffs), numberOfBatch)
    mc_configs = [MonteCarloConfiguration(div(mcConfig.Nsim, numberOfBatch), mcConfig.Nstep, mcConfig.monteCarloMethod, mcConfig.parallelMode.sub_mod) for _ = 1:numberOfBatch]
    for i = 1:numberOfBatch
        mc_configs[i] = @views @set mc_configs[i].parallelMode.seed = mcConfig.parallelMode.seeds[i]
    end
    Threads.@threads for i = 1:numberOfBatch
        price[:, i] = pricer(mcProcess, rfCurve, mc_configs[i], abstractPayoffs)
    end
    Out = sum(price, dims = 2) / numberOfBatch
    return Out
end
