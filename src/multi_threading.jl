abstract type AbstractMultiThreading <: ParallelMode end
struct MultiThreading <: AbstractMultiThreading 
numberOfThreads::Int64
MultiThreading() = new(Int64(Threads.nthreads()));
end

function pricer_macro_multithreading(model_type)
	@eval begin
		function pricer(mcProcess::$model_type,spotData::ZeroRateCurve,mcConfig::MonteCarloConfiguration{<: Integer , <: Integer , <: AbstractMonteCarloMethod ,  <: MultiThreading},abstractPayoff::AbstractPayoff)
			price=zeros(Threads.nthreads());
			Threads.@threads for i in 1:Threads.nthreads()
				price[i]=pricer(mcProcess,spotData,MonteCarloConfiguration(div(mcConfig.Nsim,Threads.nthreads()),mcConfig.Nstep,mcConfig.monteCarloMethod,SerialMode(),mcConfig.seed),abstractPayoff);
			end
			Out=sum(price)/Threads.nthreads();
			return Out;
		end
	end
end
pricer_macro_multithreading(BaseProcess)