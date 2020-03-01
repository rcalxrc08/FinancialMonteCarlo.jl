abstract type AbstractMultiThreading <: ParallelMode end
struct MultiThreading{ abstractMode <: BaseMode} <: AbstractMultiThreading 
	numberOfThreads::Int64
	sub_mod::abstractMode
	MultiThreading(sub_mod::abstractMode=SerialMode()) where { abstractMode <: BaseMode} = new{abstractMode}(Int64(Threads.nthreads()),sub_mod);
end

function pricer_macro_multithreading(model_type,payoff_type)
	@eval begin
		function pricer(mcProcess::$model_type,spotData::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration{<: Integer , <: Integer , <: AbstractMonteCarloMethod ,  <: AbstractMultiThreading , <: Random.AbstractRNG},abstractPayoff::$payoff_type)
			zero_typed=predict_output_type_zero(mcProcess,spotData,mcConfig,abstractPayoff);
			price=zeros(typeof(zero_typed),Threads.nthreads());
			Threads.@threads for i in 1:Threads.nthreads()
				price[i]=pricer(mcProcess,spotData,MonteCarloConfiguration(div(mcConfig.Nsim,Threads.nthreads()),mcConfig.Nstep,mcConfig.monteCarloMethod,mcConfig.parallelMode.sub_mod,mcConfig.seed),abstractPayoff);
			end
			Out=sum(price)/Threads.nthreads();
			return Out;
		end
	end
end

pricer_macro_multithreading(BaseProcess,AbstractPayoff)
pricer_macro_multithreading(BaseProcess,Dict{FinancialMonteCarlo.AbstractPayoff,Number})
pricer_macro_multithreading(Dict{String,FinancialMonteCarlo.AbstractMonteCarloProcess},Dict{String,Dict{FinancialMonteCarlo.AbstractPayoff,Number}})
pricer_macro_multithreading(VectorialMonteCarloProcess,Array{Dict{FinancialMonteCarlo.AbstractPayoff,Number}})

function pricer(mcProcess::BaseProcess,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration{<: Integer , <: Integer , <: AbstractMonteCarloMethod ,  <: AbstractMultiThreading, <: Random.AbstractRNG},abstractPayoffs::Array{abstractPayoff_}) where {abstractPayoff_ <: AbstractPayoff}
	zero_typed=predict_output_type_zero(mcProcess,spotData,mcConfig);
	zero_typed_vec=sum(predict_output_type_zero(abstractPayoff) for abstractPayoff in abstractPayoffs)
	price=zeros(typeof(zero_typed+zero_typed_vec),length(abstractPayoffs),Threads.nthreads());#Sigh
	Threads.@threads for i in 1:Threads.nthreads()
		price[:,i]=pricer(mcProcess,spotData,MonteCarloConfiguration(div(mcConfig.Nsim,Threads.nthreads()),mcConfig.Nstep,mcConfig.monteCarloMethod,mcConfig.parallelMode.sub_mod,mcConfig.seed),abstractPayoffs);
	end
	Out=sum(price,dims=2)/Threads.nthreads();
	return Out;
end
