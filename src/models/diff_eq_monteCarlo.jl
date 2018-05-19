
function simulate(mcProcess::MonteCarloProblem,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,T::Float64,monteCarloMode::MonteCarloMode=standard)
	
	S0=spotData.S0;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	Dt=T/Nstep
	X0=mcProcess.prob.u0[1];
	sol = solve(mcProcess,num_monte=Nsim,parallel_type=:threads,dt=Dt,adaptive=false)
	if isapprox(X0,0.0)
		S=[S0*exp.(path.u[j]) for path in sol.u, j in 1:Nstep];
		
		return S;
	else
		S=[path.u[j] for path in sol.u, j in 1:Nstep];
		
		return S;
	end
end
