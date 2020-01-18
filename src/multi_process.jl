abstract type AbstractMultiProcess <: ParallelMode end
struct MultiProcess{ abstractMode <: BaseMode}  <: AbstractMultiProcess 
	numberOfProcesses::Int64
	sub_mod::abstractMode
	MultiProcess(sub_mod::abstractMode=SerialMode()) where { abstractMode <: BaseMode} = new{abstractMode}(Int64(nworkers()),sub_mod);
end

function pricer_macro_multiprocesses(model_type,payoff_type)
	@eval begin
		function pricer(mcProcess::$model_type,spotData::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration{<: Integer , <: Integer , <: AbstractMonteCarloMethod ,  <: MultiProcess, <: Random.AbstractRNG},abstractPayoff::$payoff_type)
			#price::Number = @distributed (+) for i in 1:10
			#	Number(pricer(mcProcess,spotData,MonteCarloConfiguration(div(mcConfig.Nsim,10),mcConfig.Nstep,mcConfig.monteCarloMethod,mcConfig.parallelMode.sub_mod,mcConfig.seed+1+i),abstractPayoff));
			#end
			f(i)=pricer(mcProcess,spotData,MonteCarloConfiguration(div(mcConfig.Nsim,nworkers()),mcConfig.Nstep,mcConfig.monteCarloMethod,mcConfig.parallelMode.sub_mod,mcConfig.seed+1+i),abstractPayoff);
			Price=pmap(f,1:nworkers());
			Out=sum(Price)/nworkers();
			return Out;
		end
	end
end
pricer_macro_multiprocesses(BaseProcess,AbstractPayoff)
pricer_macro_multiprocesses(BaseProcess,Dict{FinancialMonteCarlo.AbstractPayoff,Number})
pricer_macro_multiprocesses(Dict{String,FinancialMonteCarlo.AbstractMonteCarloProcess},Dict{String,Dict{FinancialMonteCarlo.AbstractPayoff,Number}})
pricer_macro_multiprocesses(VectorialMonteCarloProcess,Array{Dict{FinancialMonteCarlo.AbstractPayoff,Number}})

function pricer(mcProcess::BaseProcess,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration{<: Integer , <: Integer , <: AbstractMonteCarloMethod ,  <: MultiProcess, <: Random.AbstractRNG},abstractPayoffs::Array{abstractPayoff_}) where {abstractPayoff_ <: AbstractPayoff}
	price = @distributed (+) for i in 1:nworkers()
		pricer(mcProcess,spotData,MonteCarloConfiguration(div(mcConfig.Nsim,nworkers()),mcConfig.Nstep,mcConfig.monteCarloMethod,mcConfig.parallelMode.sub_mod,mcConfig.seed+1+i),abstractPayoff);
	end
	Out=price./nworkers();
	return Out;
end
