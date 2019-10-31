function get_initial_u0(mcProcess::MonteCarloProblem)
	if(typeof(mcProcess.prob)<:JumpProblem)
		return  mcProcess.prob.prob.u0;
	else
		return monte_prob.prob.u0;
	end
end



mutable struct MonteCarloDiffeEqModel{num <: Number} <: ItoProcess
	model::num
	underlying::Underlying
	final_trasform::Function
	function MonteCarloDiffeEqModel(model::num,final_trasform::Function,underlying::Underlying) where {num <: DiffEqBase.AbstractEnsembleProblem}
        return new{num}(model,final_trasform,underlying)
    end
	function MonteCarloDiffeEqModel(model::num,final_trasform::Function,S0::num2) where {num <: Number,num1 <: Number, num2 <: Number}
        if σ <= 0.0
            error("Volatility must be positive")
        else
            return new{num}(model,final_trasform,Underlying(S0))
        end
    end
end


function simulate(mcProcess::MonteCarloDiffeEqModel,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	Dt=T/Nstep
	tt=collect(0.0:Dt:T)
	if T<=0.0
		error("Final time must be positive");
	end
	diffeqmodel=mcProcess.model;
	sol = solve(diffeqmodel,SOSRI(),trajectories=Nsim,parallel_type=:threads,dt=Dt,adaptive=false)
	if(!(typeof(diffeqmodel)<:JumpProblem))
		X=[path.u[j] for path in sol.u, j in 1:(Nstep+1)];
		if(get_initial_u0(diffeqmodel)==mcProcess.underlying.S0)
			return X;
		else
			return final_trasform.(X);
		end
	else
		X=[path(tt[j]) for path in sol.u, j in 1:(Nstep+1)];
		if(get_initial_u0(diffeqmodel)==mcProcess.underlying.S0)
			return X;
		else
			return final_trasform.(X);
		end
	end
end
