using .Distributed

abstract type AbstractMultiProcess <: ParallelMode end
struct MultiProcess{abstractMode <: BaseMode} <: AbstractMultiProcess
    numberOfProcesses::Int64
    nbatches::Int64
    sub_mod::abstractMode
    MultiProcess(sub_mod::abstractMode = SerialMode(), nbatches::Int64 = nworkers()) where {abstractMode <: BaseMode} = new{abstractMode}(Int64(nworkers()), nbatches, sub_mod)
end

function GeneralMonteCarloConfiguration(Nsim::num1, Nstep::num2, monteCarloMethod::abstractMonteCarloMethod, parallelMethod::baseMode, seed::Number, rng::rngType_, offset::Number) where {num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod, baseMode <: MultiProcess, rngType_ <: MersenneTwister}
    @assert Nsim > zero(num1) "Number of Simulations must be positive"
    @assert Nstep > zero(num2) "Number of Steps must be positive"
    return new{num1, num2, abstractMonteCarloMethod, baseMode, rngType_}(Nsim, Nstep, monteCarloMethod, parallelMethod, Int64(seed), div(Nsim, 2) * Nstep * (myid() == 1 ? 0 : (myid() - 2)), rng)
end

function pricer_macro_multiprocesses(model_type, payoff_type)
    @eval begin
        function pricer(mcProcess::$model_type, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:MultiProcess, <:Random.AbstractRNG}, abstractPayoff::$payoff_type)
            zero_typed = predict_output_type_zero(mcProcess, rfCurve, mcConfig, abstractPayoff)
            nbatches = mcConfig.parallelMode.nbatches
            price::typeof(zero_typed) = @sync @distributed (+) for i = 1:nbatches
                pricer(mcProcess, rfCurve, MonteCarloConfiguration(div(mcConfig.Nsim, nbatches), mcConfig.Nstep, mcConfig.monteCarloMethod, mcConfig.parallelMode.sub_mod, mcConfig.seed + 1 + i), abstractPayoff)::typeof(zero_typed)
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

function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:MultiProcess, <:Random.AbstractRNG}, abstractPayoffs::Array{abstractPayoff_}) where {abstractPayoff_ <: AbstractPayoff}
    nbatches = mcConfig.parallelMode.nbatches
    price = @distributed (+) for i = 1:nbatches
        pricer(mcProcess, rfCurve, MonteCarloConfiguration(div(mcConfig.Nsim, nbatches), mcConfig.Nstep, mcConfig.monteCarloMethod, mcConfig.parallelMode.sub_mod, mcConfig.seed + 1 + i), abstractPayoffs)
    end
    Out = price ./ nbatches
    return Out
end
