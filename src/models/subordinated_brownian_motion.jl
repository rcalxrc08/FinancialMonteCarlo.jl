type SubordinatedBrownianMotion<:AbstractMonteCarloProcess end

export SubordinatedBrownianMotion;

function simulate(mcProcess::SubordinatedBrownianMotion,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64,dt_s::Array{Float64,2},monteCarloMode::MonteCarloMode=standard)
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	if(length(mcBaseData.param)!=2)
		error("Brownian Subordinator needs 2 parameters")
	elseif(size(dt_s)!=(Nsim,Nstep))
		error("Inconsistent Time Matrix")
	end
	drift=mcBaseData.param["drift"];
	sigma=mcBaseData.param["sigma"];
	if sigma<=0.0
		error("Subordinator volatility must be positive")
	elseif T<=0.0
		error("Final time must be positive");
	end
	isDualZero=drift*sigma*dt_s[1,1]*0.0;
	X=Matrix{typeof(isDualZero)}(Nsim,Nstep+1);
	X[:,1]=isDualZero;
	if monteCarloMode==antithetic
		for i=1:Nstep
			NsimAnti=Int(floor(Nsim/2))
			Z=randn(NsimAnti);
			Z=[Z;-Z];
			# SUBORDINATED BROWNIAN MOTION (dt_s=time change)
			X[:,i+1]=X[:,i].+drift.*dt_s[:,i].+sigma.*sqrt.(dt_s[:,i]).*Z;
		end
	else
		Z=Array{Float64}(Nsim)
		for i=1:Nstep
			randn!(Z)
			# SUBORDINATED BROWNIAN MOTION (dt_s=time change)
			X[:,i+1]=X[:,i].+drift.*dt_s[:,i].+sigma.*sqrt.(dt_s[:,i]).*Z;
		end
	end
	return X;
end
