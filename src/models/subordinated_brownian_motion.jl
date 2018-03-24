type SubordinatedBrownianMotion<:AbstractMonteCarloProcess end

export SubordinatedBrownianMotion;

function simulate(mcProcess::SubordinatedBrownianMotion,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64,dt_s::Array{Float64,2})
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	if(length(mcBaseData.param)!=2)
		error("Brownian Subordinator needs 2 parameters")
	end
	drift=mcBaseData.param["drift"];
	sigma=mcBaseData.param["sigma"];
	
	isDualZero=drift*sigma*dt_s[1,1];
	X=zeros(typeof(isDualZero),Nsim,Nstep+1);
	for i=1:Nstep
		# SUBORDINATED BROWNIAN MOTION (dt_s=time change)
		X[:,i+1]=X[:,i]+drift*dt_s[:,i]+sigma*sqrt.(dt_s[:,i]).*randn(Nsim);
	end
	
	return X;
end
