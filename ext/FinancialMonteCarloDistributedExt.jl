module FinancialMonteCarloDistributedExt

using Accessors
using Distributed
using FinancialMonteCarlo: AbstractMonteCarloMethod, AbstractMonteCarloProcess, AbstractPayoff, AbstractZeroRateCurve, BaseProcess, MonteCarloConfiguration, MultiProcess, VectorialMonteCarloProcess, predict_output_type_zero
import FinancialMonteCarlo: default_nbatches, pricer

default_nbatches() = Int64(nworkers())

function pricer_macro_multiprocesses(model_type, payoff_type)
    @eval begin
        function pricer(mcProcess::$model_type, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:MultiProcess}, abstractPayoff::$payoff_type)
            zero_typed = predict_output_type_zero(mcProcess, rfCurve, mcConfig, abstractPayoff)
            nbatches = mcConfig.parallelMode.nbatches
            mc_configs = [MonteCarloConfiguration(div(mcConfig.Nsim, nbatches), mcConfig.Nstep, mcConfig.monteCarloMethod, mcConfig.parallelMode.sub_mod) for _ = 1:nbatches]
            for i = 1:nbatches
                mc_config_i = mc_configs[i]
                mc_configs[i] = @set mc_config_i.parallelMode.seed = mcConfig.parallelMode.seeds[i]
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
        mc_config_i = mc_configs[i]
        mc_configs[i] = @set mc_config_i.parallelMode.seed = mcConfig.parallelMode.seeds[i]
    end
    price = @distributed (+) for i = 1:nbatches
        pricer(mcProcess, rfCurve, mc_configs[i], abstractPayoffs)
    end
    Out = price ./ nbatches
    return Out
end

end
