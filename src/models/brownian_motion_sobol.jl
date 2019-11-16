using Sobol, StatsFuns


function simulate(mcProcess::BrownianMotion,spotData::ZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: SobolMode}
	@show "sobol mode"
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
	seq=SobolSeq(1);
	#skip(seq,Nstep*Nsim)
	#seq=SobolSeq(1);
	@inbounds for j=1:Nstep
		#seq=SobolSeq(1);
		#(j>1) ? skip(seq,j*Nsim) : nothing
		@inbounds for i=1:Nsim
			X[i,j+1]=X[i,j]+mean_bm+stddev_bm*norminvcdf(next!(seq)[1]);
		end
	end
	#@inbounds for j=1:Nstep
	#	X[:,j+1]=X[:,j].+mean_bm.+stddev_bm*norminvcdf.(next!(seq));
	#end

	return X;

end