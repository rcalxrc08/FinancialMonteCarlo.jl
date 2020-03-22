
function simulate!(X,mcProcess::HestonProcess,rfCurve::ZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode,type5},T::numb) where {numb <: Number, type1 <: Integer, type2<: Integer, type3 <: StandardMC, type5 <: Random.AbstractRNG}
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
	zero_drift=r_d(dt*0.0,dt);
	isDualZero=S0*T*σ_zero*κ*θ*λ1*σ*ρ*0.0*zero_drift;
	#X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	X[:,1].=isDualZero;
	v_m=[σ_zero+isDualZero for _ in 1:Nsim];
	for j in 1:Nstep
		e1=randn(mcBaseData.rng,Nsim);
		e2=randn(mcBaseData.rng,Nsim);
		e2=e1*ρ+e2*sqrt(1-ρ*ρ);
		tmp_=r_d((j-1)*dt,dt);
		X[:,j+1]=X[:,j]+(tmp_.-0.5*max.(v_m,0)*dt)+sqrt.(max.(v_m,0))*sqrt(dt).*e1;
		v_m+=κ_s.*(θ_s.-max.(v_m,0)).*dt+σ.*sqrt.(max.(v_m,0)).*sqrt(dt).*e2;
	end
	## Conclude
	X.=S0.*exp.(X);
	return ;

end

function simulate!(X,mcProcess::HestonProcess,rfCurve::ZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,SerialMode,type5},T::numb) where {numb <: Number, type1 <: Integer, type2<: Integer, type3 <: AntitheticMC, type5 <: Random.AbstractRNG}
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

	dt=T/Nstep
	isDualZero=S0*T*σ_zero*κ*θ*λ1*σ*ρ*0.0;
	#X=Matrix{typeof(isDualZero)}(undef,Nsim,Nstep+1);
	X[:,1].=isDualZero;
	v_m=[σ_zero+isDualZero for _ in 1:Nsim];
	Nsim_2=div(Nsim,2)
	r_d=r-d;

	for j in 1:Nstep
		e1=randn(mcBaseData.rng,Nsim_2);
		e2=randn(mcBaseData.rng,Nsim_2);
		e1=[e1;-e1];
		e2=[e2;-e2];
		e2=e1.*ρ.+e2.*sqrt(1-ρ*ρ);
		tmp_=r_d((j-1)*dt,dt);
		X[:,j+1]=X[:,j]+(tmp_.-0.5.*max.(v_m,0).*dt)+sqrt.(max.(v_m,0)).*sqrt(dt).*e1;
		v_m+=κ_s.*(θ_s.-max.(v_m,0)).*dt+(σ*sqrt(dt)).*sqrt.(max.(v_m,0)).*e2;
	end
	
	## Conclude
	X.=S0.*exp.(X);
	return;

end