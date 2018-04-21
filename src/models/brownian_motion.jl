type BrownianMotion<:ItoProcess end

export BrownianMotion;

function simulate(mcProcess::BrownianMotion,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64)
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	if(length(mcBaseData.param)!=2)
		error("Brownian Motion needs 2 parameters")
	end
	volatility=mcBaseData.param["sigma"];
	drift=mcBaseData.param["drift"];
	dt=T/Nstep
	meanW=drift*dt
	varW=volatility*sqrt(dt)
	isDualZero=meanW*varW*0.0;
	X=Matrix{typeof(isDualZero)}(Nsim,Nstep+1);
	X[:,1]=isDualZero;
	for j in 1:Nstep
		X[:,j+1]=X[:,j]+meanW.+varW.*randn(Nsim);
	end

	return X;

end
