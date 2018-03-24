type NormalInverseGaussianProcess<:AbstractMonteCarloProcess end

export NormalInverseGaussianProcess;

function simulate(mcProcess::NormalInverseGaussianProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64)
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
	
	psi1(v::Number)::Number=(1-sqrt(1.0+ ((sigma*v)^2-2.0*1im*theta1*v)*k1))/k1;
	drift=r-d-real(psi1(-1im));
	const dict1=Dict{String,Number}("sigma"=>sigma, "drift" => drift)
	subParam=MonteCarloBaseData(dict1,Nsim,Nstep);
	
	#Simulate subordinator
	IGRandomVariable=InverseGaussian(dt,dt*dt/k1);
	dt_s=quantile.(IGRandomVariable,rand(Nsim,Nstep));
	
	X=simulate(SubordinatedBrownianMotion(),spotData,subParam,T,dt_s);

	S=S0*exp.(X);
	
	return S;
	
end
