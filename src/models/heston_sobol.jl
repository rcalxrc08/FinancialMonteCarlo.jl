using Sobol, StatsFuns

function simulate!(X,mcProcess::HestonProcess,rfCurve::ZeroRate,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode,type5},T::numb) where {numb <: Number, type1 <: Integer, type2<: Integer, type3 <: SobolMode, type5 <: Random.AbstractRNG}
	r=rfCurve.r;
	S0=mcProcess.underlying.S0;
	d=dividend(mcProcess);
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ=mcProcess.σ;
	σ_zero=mcProcess.σ_zero;
	λ1=mcProcess.λ;
	κ=mcProcess.κ;
	ρ=mcProcess.ρ;
	θ=mcProcess.θ;
	@assert T>0.0

	####Simulation
	## Simulate
	κ_s=κ+λ1;
	θ_s=κ*θ/κ_s;
	seq=SobolSeq(2*Nstep);
	skip(seq,Nstep*Nsim)
	dt=T/Nstep
	isDualZero=T*r*σ_zero*θ_s*κ_s*σ*ρ*0.0*S0;
	view(X,:,1).=isDualZero;
	for i in 1:Nsim
		v_m=σ_zero^2;
		vec=norminvcdf.(next!(seq))
		for j in 1:Nstep
			@views e1=vec[j];
			@views e2=e1*ρ+vec[Nstep+j]*sqrt(1-ρ*ρ);
			@views X[i,j+1]=X[i,j]+((r-d)-0.5*v_m)*dt+sqrt(v_m)*sqrt(dt)*e1;
			v_m+=κ_s*(θ_s-v_m)*dt+σ*sqrt(v_m)*sqrt(dt)*e2;
			#when v_m = 0.0, the derivative becomes NaN
			v_m=max(v_m,isDualZero);
		end
	end
	## Conclude
	X.=S0.*exp.(X);
	return;

end



function simulate!(X,mcProcess::HestonProcess,rfCurve::ZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode,type5},T::numb) where {numb <: Number, type1 <: Integer, type2<: Integer, type3 <: SobolMode, type5 <: Random.AbstractRNG}
	r=rfCurve.r;
	S0=mcProcess.underlying.S0;
	d=dividend(mcProcess);
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ=mcProcess.σ;
	σ_zero=mcProcess.σ_zero;
	λ1=mcProcess.λ;
	κ=mcProcess.κ;
	ρ=mcProcess.ρ;
	θ=mcProcess.θ;
	@assert T>0.0

	####Simulation
	## Simulate
	κ_s=κ+λ1;
	θ_s=κ*θ/(κ+λ1);
	r_d=r-d;
	dt=T/Nstep
	seq=SobolSeq(2*Nstep);
	skip(seq,Nstep*Nsim)
	zero_drift=r_d(dt*0.0,dt);
	isDualZero=S0*T*σ_zero*κ*θ*λ1*σ*ρ*0.0*zero_drift;
	#X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	X[:,1].=isDualZero;
	tmp_=[r_d((j-1)*dt,dt) for j in 1:Nstep];
	for i in 1:Nsim
		v_m=σ_zero^2;
		vec=norminvcdf.(next!(seq))
		
		for j in 1:Nstep
			@views e1=vec[j];
			@views e2=e1*ρ+vec[Nstep+j]*sqrt(1-ρ*ρ);
			@views X[i,j+1]=X[i,j]+(tmp_[j]-0.5*v_m)*dt+sqrt(v_m)*sqrt(dt)*e1;
			v_m+=κ_s*(θ_s-v_m)*dt+σ*sqrt(v_m)*sqrt(dt)*e2;
			#when v_m = 0.0, the derivative becomes NaN
			v_m=max(v_m,isDualZero);
		end
	end
	## Conclude
	X.=S0.*exp.(X);
	return ;

end