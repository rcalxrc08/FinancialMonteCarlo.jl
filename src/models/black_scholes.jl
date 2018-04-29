type BlackScholesProcess<:ItoProcess end

export BlackScholesProcess;

function simulate(mcProcess::BlackScholesProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64,mode1::MonteCarloMode=standard)
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	if(length(mcBaseData.param)!=1)
		error("Black Scholes Model needs 1 parameters")
	end
	sigma_gbm=mcBaseData.param["sigma"];
	mu_gbm=r-d;
	const dictGBM=Dict{String,Number}("sigma"=>sigma_gbm, "drift" => mu_gbm)
	GeomData=MonteCarloBaseData(dictGBM,mcBaseData.Nsim,mcBaseData.Nstep)
	
	S=S0.*simulate(GeometricBrownianMotion(),spotData,GeomData,T,mode1)
	
	return S;
	
end