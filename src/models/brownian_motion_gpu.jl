
function simulate(mcProcess::BrownianMotion,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::numb,monteCarloMode::MonteCarloMode,parallelMode::CudaMode) where {numb<:Number}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;

	σ=mcProcess.σ;
	μ=mcProcess.μ;
	if T<=0.0
		error("Final time must be positive");
	end
	dt=T/Nstep
	mean_bm=μ*dt
	stddev_bm=σ*sqrt(dt)
	if monteCarloMode==antithetic
		Nsim_2=Int(floor(Nsim/2))
		if Nsim_2*2!=Nsim
			error("Antithetic support only odd number of simulations")
		end
		mean_bm_f=Float32(mean_bm);
		stddev_bm_f=Float32(stddev_bm);
		isDualZero=mean_bm_f*stddev_bm_f*zero(Float32);
		X_cu=CuMatrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
		for i=1:Nstep
			Z=randn(Nsim_2);
			Z=[Z;-Z];
			X_cu[:,i+1]=X_cu[:,i]+(mean_bm_f.+stddev_bm_f.*cu(Z));
		end
		X=Matrix(X_cu);
		return X;
	else
		mean_bm_f=Float32(mean_bm);
		stddev_bm_f=Float32(stddev_bm);
		isDualZero=mean_bm_f*stddev_bm_f*zero(Float32);
		X_cu=CuMatrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
		for i=1:Nstep
			X_cu[:,i+1]=X_cu[:,i]+(mean_bm_f.+stddev_bm_f.*cu(randn(Float32,Nsim)));
		end
		X=Matrix(X_cu);
		return X;
	end

end
