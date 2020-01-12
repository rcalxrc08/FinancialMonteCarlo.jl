abstract type AbstractMultiThreading <: ParallelMode end
struct MultiThreading <: AbstractMultiThreading 
	numberOfThreads::Int64
	MultiThreading() = new(Int64(Threads.nthreads()));
end

function pricer_macro_multithreading(model_type,payoff_type)
	@eval begin
		function pricer(mcProcess::$model_type,spotData::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration{<: Integer , <: Integer , <: AbstractMonteCarloMethod ,  <: MultiThreading},abstractPayoff::$payoff_type)
			price=zeros(Number,Threads.nthreads());#Sigh
			Threads.@threads for i in 1:Threads.nthreads()
				price[i]=pricer(mcProcess,spotData,MonteCarloConfiguration(div(mcConfig.Nsim,Threads.nthreads()),mcConfig.Nstep,mcConfig.monteCarloMethod,SerialMode(),mcConfig.seed),abstractPayoff);
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

function pricer(mcProcess::BaseProcess,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration{<: Integer , <: Integer , <: AbstractMonteCarloMethod ,  <: MultiThreading},abstractPayoffs::Array{abstractPayoff_}) where {abstractPayoff_ <: AbstractPayoff}
	price=zeros(Number,length(abstractPayoffs),Threads.nthreads());#Sigh
	Threads.@threads for i in 1:Threads.nthreads()
		price[:,i]=pricer(mcProcess,spotData,MonteCarloConfiguration(div(mcConfig.Nsim,Threads.nthreads()),mcConfig.Nstep,mcConfig.monteCarloMethod,SerialMode(),mcConfig.seed),abstractPayoffs);
	end
	Out=sum(price,dims=2)/Threads.nthreads();
	return Out;
end
