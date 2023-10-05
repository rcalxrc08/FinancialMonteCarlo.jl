using .Distributed

abstract type AbstractMultiProcess <: ParallelMode end
struct MultiProcess{abstractMode <: BaseMode} <: AbstractMultiProcess
    nbatches::Int64
    sub_mod::abstractMode
    seeds::Array{Int64}
    function MultiProcess(sub_mod::abstractMode = SerialMode(), nbatches::Int64 = Int64(nworkers()), seeds::Array{Int64} = collect(1:nbatches)) where {abstractMode <: BaseMode}
        @assert length(seeds) == nbatches
        return new{abstractMode}(nbatches, sub_mod, seeds)
    end
end

function pricer_macro_multiprocesses(model_type, payoff_type)
    @eval begin
        function pricer(mcProcess::$model_type, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:MultiProcess}, abstractPayoff::$payoff_type)
            zero_typed = predict_output_type_zero(mcProcess, rfCurve, mcConfig, abstractPayoff)
            nbatches = mcConfig.parallelMode.nbatches
            mc_configs = [MonteCarloConfiguration(div(mcConfig.Nsim, nbatches), mcConfig.Nstep, mcConfig.monteCarloMethod, mcConfig.parallelMode.sub_mod) for _ = 1:nbatches]
            for i = 1:nbatches
                @set mc_configs[i].parallelMode.seed = mcConfig.parallelMode.seeds[i]
            end
            price::typeof(zero_typed) = @sync @distributed (+) for i = 1:nbatches
                pricer(mcProcess, rfCurve, mc_configs[i], abstractPayoff)::typeof(zero_typed)
            end
            Out = price / nbatches
            return Out
        end
    end
end
pricer_macro_multiprocesses(BaseProcess, AbstractPayoff)
pricer_macro_multiprocesses(BaseProcess, Dict{AbstractPayoff, Number})
pricer_macro_multiprocesses(Dict{String, AbstractMonteCarloProcess}, Dict{String, Dict{AbstractPayoff, Number}})
pricer_macro_multiprocesses(VectorialMonteCarloProcess, Array{Dict{AbstractPayoff, Number}})

function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:MultiProcess}, abstractPayoffs::Array{abstractPayoff_}) where {abstractPayoff_ <: AbstractPayoff}
    nbatches = mcConfig.parallelMode.nbatches
    mc_configs = [MonteCarloConfiguration(div(mcConfig.Nsim, nbatches), mcConfig.Nstep, mcConfig.monteCarloMethod, mcConfig.parallelMode.sub_mod) for _ = 1:nbatches]
    for i = 1:nbatches
        @set mc_configs[i].parallelMode.seed = mcConfig.parallelMode.seeds[i]
    end
    price = @distributed (+) for i = 1:nbatches
        pricer(mcProcess, rfCurve, mc_configs[i], abstractPayoffs)
    end
    Out = price ./ nbatches
    return Out
end
