using Sobol, StatsFuns


function simulate(mcProcess::BrownianMotion,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: SobolMode}
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
	isDualZero=mean_bm*stddev_bm*0.0+mcProcess.underlying.S0;
	X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	X[:,1].=isDualZero;
	seq=SobolSeq(Nstep);
	skip(seq,Nstep*Nsim)
	@inbounds for i=1:Nsim
		vec=norminvcdf.(next!(seq))
		@inbounds for j=1:Nstep
			@views X[i,j+1]=X[i,j]+mean_bm+stddev_bm*vec[j];
		end
	end
	return X;
end