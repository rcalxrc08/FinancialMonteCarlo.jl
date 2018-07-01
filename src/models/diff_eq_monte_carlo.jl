
function simulate(mcProcess::MonteCarloProblem,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::numb,monteCarloMode::MonteCarloMode=standard) where {numb<:Number}
	
	S0=spotData.S0;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	Dt=T/Nstep
	
	if T<=0.0
		error("Final time must be positive");
	end
	sol = solve(mcProcess,SOSRI(),num_monte=Nsim,parallel_type=:threads,dt=Dt,adaptive=false)
	X0=sol.u[1].u[1];
	if isapprox(X0,0.0)
		S=[S0*exp.(path.u[j]) for path in sol.u, j in 1:Nstep];
		
		return S;
	else
		S=[path.u[j] for path in sol.u, j in 1:Nstep];
		
		return S;
	end
end
