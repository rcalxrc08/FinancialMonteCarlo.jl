
function simulate(mcProcess::finiteActivityProcess,rfCurve::ZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where { finiteActivityProcess <: FiniteActivityProcess, numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	r=rfCurve.r;
	S0=mcProcess.underlying.S0;
	d=dividend(mcProcess);
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ=mcProcess.σ;
	λ1=mcProcess.λ;
	
	if T<=0.0
		error("Final time must be positive");
	end

	####Simulation
	## Simulate
	# r-d-psi(-i)
	#drift_RN=r-d-σ^2/2-λ1*(p/(λ₊-1)-(1-p)/(λ₋+1));
	drift_RN=r-d-compute_drift(mcProcess);
	X=Matrix(simulate(BrownianMotion(σ,drift_RN,Underlying(0.0)),rfCurve,mcBaseData,T))

	t=range(0.0, stop=T, length=Nstep+1);
	PoissonRV=Poisson(λ1*T);
	NJumps=quantile.(PoissonRV,rand(mcBaseData.rng,Nsim));

	for ii in 1:Nsim
		Njumps_=NJumps[ii];
		# Simulate the times of jump
		tjumps=sort(rand(mcBaseData.rng,Njumps_)*T);
		for tjump in tjumps
			# Add the jump size
			idx1=findfirst(x->x>=tjump,t);
			jump_size=compute_jump_size(mcProcess,mcBaseData);
			X[ii,idx1:end].+=jump_size; #add jump component
			
		end
	end
	## Conclude
	S=S0.*exp.(X);
	return S;

end