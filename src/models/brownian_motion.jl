type BrownianMotion<:ItoProcess end

export BrownianMotion;

function simulate(mcProcess::BrownianMotion,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64,monteCarloMode::MonteCarloMode=standard)
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	if(length(mcBaseData.param)!=2)
		error("Brownian Motion needs 2 parameters")
	end
	sigma=mcBaseData.param["sigma"];
	drift=mcBaseData.param["drift"];
	if sigma<=0.0
		error("Brownian motion volatility must be positive");
	elseif T<=0.0
		error("Final time must be positive");
	end
	dt=T/Nstep
	mean_bm=drift*dt
	stddev_bm=sigma*sqrt(dt)
	isDualZero=mean_bm*stddev_bm*0.0;
	X=Matrix{typeof(isDualZero)}(Nsim,Nstep+1);
	X[:,1]=isDualZero;
	if monteCarloMode==antithetic
		Nsim_2=Int(floor(Nsim/2))
		if Nsim_2*2!=Nsim
			error("Antithetic support only odd number of simulations")
		end
		for j in 1:Nstep
			Z=randn(Nsim_2);
			Z=[Z;-Z];
			X[:,j+1]=X[:,j]+mean_bm.+stddev_bm.*Z;
		end
	else
		Z=Array{Float64}(Nsim)
		for j in 1:Nstep
			randn!(Z)
			X[:,j+1]=X[:,j]+mean_bm.+stddev_bm.*Z;
		end
	end

	return X;

end
