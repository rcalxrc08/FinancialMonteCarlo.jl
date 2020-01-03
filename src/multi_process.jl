abstract type AbstractMultiProcess <: ParallelMode end
struct MultiProcess <: AbstractMultiProcess 
numberOfProcesses::Int64
MultiProcess() = new(Int64(nworkers()));
end

function pricer_macro_multiprocesses(model_type)
	@eval begin
		function pricer(mcProcess::$model_type,spotData::ZeroRate,mcConfig::MonteCarloConfiguration{<: Integer , <: Integer , <: AbstractMonteCarloMethod ,  <: MultiProcess},abstractPayoff::AbstractPayoff)
			price = @distributed (+) for i in 1:nworkers()
				pricer(mcProcess,spotData,MonteCarloConfiguration(div(mcConfig.Nsim,nworkers()),mcConfig.Nstep,mcConfig.monteCarloMethod,SerialMode(),mcConfig.seed+1+i),abstractPayoff);
			end
			Out=price/nworkers();
			return Out;
		end
	end
end
pricer_macro_multiprocesses(BaseProcess)