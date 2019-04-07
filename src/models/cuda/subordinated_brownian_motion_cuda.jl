
function simulate(mcProcess::SubordinatedBrownianMotion,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::numb,dt_s::Array{num1,2},monteCarloMode::MonteCarloMode,parallelMode::CudaMode) where {numb,num1<:Number}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	if(size(dt_s)!=(Nsim,Nstep))
		error("Inconsistent Time Matrix")
	end
	drift=mcProcess.drift;
	sigma=mcProcess.sigma;
	if T<=0.0
		error("Final time must be positive");
	end
	isDualZero=drift*sigma*dt_s[1,1]*zero(Float32);
	X_cu=CuMatrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	X[:,1].=isDualZero;
	if monteCarloMode==antithetic
		for i=1:Nstep
			NsimAnti=Int(floor(Nsim/2))
			Z=CuArrays.CURAND.curandn(Float32,NsimAnti));
			Z=[Z;-Z];
			# SUBORDINATED BROWNIAN MOTION (dt_s=time change)
			X[:,i+1]=X[:,i].+drift.*dt_s[:,i].+sigma.*sqrt.(dt_s[:,i]).*Z;
		end
	else
		Z=Array{Float64}(undef,Nsim)
		for i=1:Nstep
			# SUBORDINATED BROWNIAN MOTION (dt_s=time change)
			X[:,i+1]=X[:,i].+drift.*dt_s[:,i].+sigma.*sqrt.(dt_s[:,i]).*CuArrays.CURAND.curandn(Float32,Nsim));
		end
	end
	return X;
end
