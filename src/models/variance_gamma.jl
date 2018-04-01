type VarianceGammaProcess<:AbstractMonteCarloProcess end

export VarianceGammaProcess;

function simulate(mcProcess::VarianceGammaProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64)
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	if(length(mcBaseData.param)!=3)
		error("Variance Gamma Model needs 3 parameters")
	end
	sigma=mcBaseData.param["sigma"];
	theta1=mcBaseData.param["theta"];
	k1=mcBaseData.param["k"];
	dt=T/Nstep;
	psi1(u::Number)::Number=-1/k1*log(1+u*u*sigma*sigma*k1/2.0-1im*theta1*k1*u);
	drift=r-d-real(psi1(-1im));
	
	gammaRandomVariable=Gamma(dt/k1);
	dt_s=k1.*quantile.(gammaRandomVariable,rand(Nsim,Nstep));
	const dict1=Dict{String,Number}("sigma"=>sigma, "drift" => drift)
	subParam=MonteCarloBaseData(dict1,Nsim,Nstep);
	
	X=simulate(SubordinatedBrownianMotion(),spotData,subParam,T,dt_s);

	S=S0.*exp.(X);
	
	return S;
	
end
