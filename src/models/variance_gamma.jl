type VarianceGammaProcess<:InfiniteActivityProcess end

export VarianceGammaProcess;

function simulate(mcProcess::VarianceGammaProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64,monteCarloMode::MonteCarloMode=standard)
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	sigma=mcBaseData.param["sigma"];
	theta1=mcBaseData.param["theta"];
	k1=mcBaseData.param["k"];
	if(length(mcBaseData.param)!=3)
		error("Variance Gamma Model needs 3 parameters")
	elseif 1-sigma*sigma*k1/2.0-theta1*k1<=0.0
		error("Parameters with unfeasible values")
	elseif k1<=0.0
		error("k must be positive");
	elseif T<=0.0
		error("Final time must be positive");
	end
	
	dt=T/Nstep;
	#-1/p[3]*log(1+u*u*p[1]*p[1]*p[3]/2.0-1im*p[2]*p[3]*u);
	psi1=-1/k1*log(1-sigma*sigma*k1/2.0-theta1*k1);
	#1-sigma*sigma*k1/2.0-theta1*k1
	drift=r-d-psi1;
	
	gammaRandomVariable=Gamma(dt/k1);
	dt_s=k1.*quantile.(gammaRandomVariable,rand(Nsim,Nstep));
	const dict1=Dict{String,Number}("sigma"=>sigma, "drift" => drift)
	subParam=MonteCarloBaseData(dict1,Nsim,Nstep);
	
	X=simulate(SubordinatedBrownianMotion(),spotData,subParam,T,dt_s,monteCarloMode);

	S=S0.*exp.(X);
	
	return S;
	
end
