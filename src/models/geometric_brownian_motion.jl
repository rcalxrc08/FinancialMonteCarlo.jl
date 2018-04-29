type GeometricBrownianMotion<:ItoProcess end

function simulate(mcProcess::GeometricBrownianMotion,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64,mode1::MonteCarloMode=standard)
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	if(length(mcBaseData.param)!=2)
		error("Geometric Brownian Motion needs 2 parameters")
	end
	sigma_gbm=mcBaseData.param["sigma"];
	mu_gbm=mcBaseData.param["drift"];
	drift_bm=mu_gbm-sigma_gbm^2/2;
	const dictGBM=Dict{String,Number}("sigma"=>sigma_gbm, "drift" => drift_bm)
	BrownianData=MonteCarloBaseData(dictGBM,mcBaseData.Nsim,mcBaseData.Nstep)
	X=simulate(BrownianMotion(),spotData,BrownianData,T,mode1)
	S=exp.(X);
	return S;
end