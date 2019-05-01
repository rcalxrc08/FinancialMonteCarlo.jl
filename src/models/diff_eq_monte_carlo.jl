
function simulate(mcProcess::MonteCarloProblem,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::numb,monteCarloMode::MonteCarloMode=standard,parallelMode::BaseMode=SerialMode()) where {numb<:Number}
	
	S0=spotData.S0;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	Dt=T/Nstep
	tt=collect(0.0:Dt:T)
	if T<=0.0
		error("Final time must be positive");
	end
	sol = solve(mcProcess,SOSRI(),num_monte=Nsim,parallel_type=:threads,dt=Dt,adaptive=false)
	X0=sol.u[1].u[1];
	if(!(typeof(mcProcess.prob)<:JumpProblem))
		if isapprox(X0,0.0)
			S=[S0*exp.(path.u[j]) for path in sol.u, j in 1:(Nstep+1)];
			
			return S;
		else
			S=[path.u[j] for path in sol.u, j in 1:(Nstep+1)];
			
			return S;
		end
	else
		if isapprox(X0,0.0)
			S=[S0*exp.(path(tt[j])) for path in sol.u, j in 1:(Nstep+1)];
			
			return S;
		else
			S=[path(tt[j]) for path in sol.u, j in 1:(Nstep+1)];
			
			return S;
		end
	end
end
