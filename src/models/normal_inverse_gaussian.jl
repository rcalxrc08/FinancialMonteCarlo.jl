type NormalInverseGaussianProcess<:InfiniteActivityProcess end

export NormalInverseGaussianProcess;

function simulate(mcProcess::NormalInverseGaussianProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64,monteCarloMode::MonteCarloMode=standard)
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
	elseif 1.0-(sigma^2+2.0*theta1)*k1<0.0
		error("Parameters with unfeasible values")
	elseif k1<=0.0
		error("k must be positive");
	elseif T<=0.0
		error("Final time must be positive");
	end
	
	dt=T/Nstep;
	#psi1(v::Number)::Number=(1-sqrt(1.0+ ((sigma*v)^2-2.0*1im*theta1*v)*k1))/k1;
	psi1=(1-sqrt(1.0-(sigma^2+2.0*theta1)*k1))/k1;
	drift=r-d-psi1;
	const dict1=Dict{String,Number}("sigma"=>sigma, "drift" => drift)
	subParam=MonteCarloBaseData(dict1,Nsim,Nstep);
	
	#Simulate subordinator
	IGRandomVariable=InverseGaussian(dt,dt*dt/k1);
	dt_s=quantile.(IGRandomVariable,rand(Nsim,Nstep));
	
	X=simulate(SubordinatedBrownianMotion(),spotData,subParam,T,dt_s,monteCarloMode);

	S=S0.*exp.(X);
	
	return S;
	
end
