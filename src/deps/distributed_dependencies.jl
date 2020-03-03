using .Distributed

abstract type AbstractMultiProcess <: ParallelMode end
struct MultiProcess{ abstractMode <: BaseMode}  <: AbstractMultiProcess 
	numberOfProcesses::Int64
	ns::Int64
	sub_mod::abstractMode
	MultiProcess(sub_mod::abstractMode=SerialMode(),ns::Int64=nworkers()) where { abstractMode <: BaseMode} = new{abstractMode}(Int64(nworkers()),ns,sub_mod);
end

function GeneralMonteCarloConfiguration(Nsim::num1,Nstep::num2,monteCarloMethod::abstractMonteCarloMethod,parallelMethod::baseMode,seed::Number,rng::rngType_,offset::Number) where {num1 <: Integer, num2 <: Integer, abstractMonteCarloMethod <: AbstractMonteCarloMethod, baseMode <: MultiProcess, rngType_ <: MersenneTwister}
	if Nsim <= zero(num1)
		error("Number of Simulations must be positive")
	elseif Nstep <= zero(num2)
		error("Number of Steps must be positive")
	else
		return new{num1,num2,abstractMonteCarloMethod,baseMode,rngType_}(Nsim,Nstep,monteCarloMethod,parallelMethod,Int64(seed),div(Nsim,2)*Nstep*( myid() == 1 ? 0 : ( myid()-2) ),rng)
	end
end


function pricer_macro_multiprocesses(model_type,payoff_type)
	@eval begin
		function pricer(mcProcess::$model_type,spotData::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration{<: Integer , <: Integer , <: AbstractMonteCarloMethod ,  <: MultiProcess, <: Random.AbstractRNG},abstractPayoff::$payoff_type)
			zero_typed=predict_output_type_zero(mcProcess,spotData,mcConfig,abstractPayoff);
			ns=mcConfig.parallelMode.ns;
			#ns=100;
			price::typeof(zero_typed) =@sync @distributed (+) for i in 1:ns
				pricer(mcProcess,spotData,MonteCarloConfiguration(div(mcConfig.Nsim,ns),mcConfig.Nstep,mcConfig.monteCarloMethod,mcConfig.parallelMode.sub_mod,mcConfig.seed+1+i),abstractPayoff)::typeof(zero_typed);
			end
			Out=price/ns;
			#f(i)=pricer(mcProcess,spotData,MonteCarloConfiguration(div(mcConfig.Nsim,nworkers()),mcConfig.Nstep,mcConfig.monteCarloMethod,mcConfig.parallelMode.sub_mod,mcConfig.seed+1+i),abstractPayoff);
			#Price=pmap(f,1:nworkers());
			#Out=sum(Price)/nworkers();
			return Out;
		end
	end
end
pricer_macro_multiprocesses(BaseProcess,AbstractPayoff)
pricer_macro_multiprocesses(BaseProcess,Dict{AbstractPayoff,Number})
pricer_macro_multiprocesses(Dict{String,AbstractMonteCarloProcess},Dict{String,Dict{AbstractPayoff,Number}})
pricer_macro_multiprocesses(VectorialMonteCarloProcess,Array{Dict{AbstractPayoff,Number}})

function pricer(mcProcess::BaseProcess,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration{<: Integer , <: Integer , <: AbstractMonteCarloMethod ,  <: MultiProcess, <: Random.AbstractRNG},abstractPayoffs::Array{abstractPayoff_}) where {abstractPayoff_ <: AbstractPayoff}
	price = @distributed (+) for i in 1:nworkers()
		pricer(mcProcess,spotData,MonteCarloConfiguration(div(mcConfig.Nsim,nworkers()),mcConfig.Nstep,mcConfig.monteCarloMethod,mcConfig.parallelMode.sub_mod,mcConfig.seed+1+i),abstractPayoff);
	end
	Out=price./nworkers();
	return Out;
end
