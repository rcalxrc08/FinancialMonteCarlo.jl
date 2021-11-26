
function simulate!(X,mcProcess::finiteActivityProcess,rfCurve::AbstractZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4,type5},T::numb) where { finiteActivityProcess <: FiniteActivityProcess, numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode, type5 <: Random.AbstractRNG}
	r=rfCurve.r;
	d=dividend(mcProcess);
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ=mcProcess.σ;
	λ1=mcProcess.λ;
	
	@assert T>0.0

	####Simulation
	## Simulate
	# r-d-psi(-i)
	drift_RN=(r-d)-σ^2/2-λ*compute_drift(mcProcess);
	simulate!(X,BrownianMotion(σ,drift_RN),mcBaseData,T)

	dt=T/Nstep;
	PoissonRV=Poisson(λ1*T);
	NJumps=rand(mcBaseData.rng,PoissonRV,Nsim);

	for (index_, njumps_) in enumerate(NJumps)
		# Simulate the times of jump
		@views tjumps=rand(mcBaseData.rng,njumps_)*T;
		for tjump in tjumps
			# Add the jump size
			idx1=ceil(Int64,tjump/dt)+1;
			#jump_size=compute_jump_size(mcProcess,mcBaseData);
			@views X[index_,idx1:end].+=compute_jump_size(mcProcess,mcBaseData); #add jump component
			
		end
	end
	## Conclude
	S0=mcProcess.underlying.S0;
	#f(x)=S0*exp(x);
	X.=S0.*exp.(X)
	#broadcast!(f,X,X)
	nothing;

end
