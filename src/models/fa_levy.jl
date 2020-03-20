
function simulate!(X,mcProcess::finiteActivityProcess,rfCurve::AbstractZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4,type5},T::numb) where { finiteActivityProcess <: FiniteActivityProcess, numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode, type5 <: Random.AbstractRNG}
	r=rfCurve.r;
	S0=mcProcess.underlying.S0;
	d=dividend(mcProcess);
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ=mcProcess.σ;
	λ1=mcProcess.λ;
	
	@assert T>0.0

	####Simulation
	## Simulate
	# r-d-psi(-i)
	drift_RN=(r-d)-compute_drift(mcProcess);
	simulate!(X,BrownianMotion(σ,drift_RN,Underlying(0.0)),mcBaseData,T)

	t=collect(range(zero(T), stop=T, length=Nstep+1));
	PoissonRV=Poisson(λ1*T);
	NJumps=rand(mcBaseData.rng,PoissonRV,Nsim);

	for ii in eachindex(NJumps)
		# Simulate the times of jump
		@views tjumps=rand(mcBaseData.rng,NJumps[ii])*T;
		for tjump in tjumps
			# Add the jump size
			idx1=findfirst(x->x>=tjump,t);
			#jump_size=compute_jump_size(mcProcess,mcBaseData);
			@views X[ii,idx1:end].+=compute_jump_size(mcProcess,mcBaseData); #add jump component
			
		end
	end
	## Conclude
	S0=mcProcess.underlying.S0;
	#f(x)=S0*exp(x);
	X.=S0.*exp.(X)
	#broadcast!(f,X,X)
	nothing;

end
