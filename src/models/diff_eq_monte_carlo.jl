
function simulate(mcProcess::MonteCarloProblem,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	Dt=T/Nstep
	tt=collect(0.0:Dt:T)
	if T<=0.0
		error("Final time must be positive");
	end
	sol = solve(mcProcess,SOSRI(),trajectories=Nsim,parallel_type=:threads,dt=Dt,adaptive=false)
	X0=sol.u[1].u[1];
	if(!(typeof(mcProcess.prob)<:JumpProblem))
		S=[path.u[j] for path in sol.u, j in 1:(Nstep+1)];
		return S;
	else
		S=[path(tt[j]) for path in sol.u, j in 1:(Nstep+1)];
		return S;
	end
end
