
function simulate(mcProcess::finiteActivityProcess,rfCurve::AbstractZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4,type5},T::numb) where { finiteActivityProcess <: FiniteActivityProcess, numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode, type5 <: Random.AbstractRNG}
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
	X=Matrix(simulate(BrownianMotion(σ,drift_RN,Underlying(0.0)),mcBaseData,T))

	t=collect(range(0.0, stop=T, length=Nstep+1));
	PoissonRV=Poisson(λ1*T);
	NJumps=quantile.(PoissonRV,rand(mcBaseData.rng,Nsim));

	for ii in 1:Nsim
		# Simulate the times of jump
		@views tjumps=rand(mcBaseData.rng,NJumps[ii])*T;
		for tjump in tjumps
			# Add the jump size
			idx1=findfirst(x->x>=tjump,t);
			@views X[ii,idx1:end].+=compute_jump_size(mcProcess,mcBaseData); #add jump component
			
		end
	end
	## Conclude
	S=S0.*exp.(X);
	return S;

end
