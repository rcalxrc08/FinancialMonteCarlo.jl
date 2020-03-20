using Sobol, StatsFuns


function simulate!(X,mcProcess::BrownianMotion,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode,type5},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: SobolMode, type5 <: Random.AbstractRNG}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ=mcProcess.σ;
	μ=mcProcess.μ;
	@assert T>0.0
	dt=T/Nstep
	mean_bm=μ*dt
	stddev_bm=σ*sqrt(dt)
	isDualZero=mean_bm*stddev_bm*0.0+mcProcess.underlying.S0;
	view(X,:,1).=isDualZero;
	seq=SobolSeq(Nstep);
	skip(seq,Nstep*Nsim)
	@inbounds for i=1:Nsim
		vec=norminvcdf.(next!(seq))
		@inbounds for j=1:Nstep
			@views X[i,j+1]=X[i,j]+mean_bm+stddev_bm*vec[j];
		end
	end
	nothing
end


function simulate!(X,mcProcess::BrownianMotionVec,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode,type5},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: SobolMode, type5 <: Random.AbstractRNG}
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ=mcProcess.σ;
	μ=mcProcess.μ;
	@assert T>0.0
	dt=T/Nstep
	stddev_bm=σ*sqrt(dt)
	zero_drift=μ(dt*0.0,dt);
	isDualZero=stddev_bm*0.0*mcProcess.underlying.S0*zero_drift;
	view(X,:,1).=isDualZero;	
	seq=SobolSeq(Nstep);
	skip(seq,Nstep*Nsim)
	tmp_=[μ((j-1)*dt,dt) for j in 1:Nstep];
	@inbounds for i=1:Nsim
		vec=norminvcdf.(next!(seq))
		@inbounds for j=1:Nstep
			@views X[i,j+1]=X[i,j]+tmp_[j]+stddev_bm*vec[j];
		end
	end

	nothing

end