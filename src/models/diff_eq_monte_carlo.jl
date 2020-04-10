mutable struct MonteCarloDiffEqModel{absdiffeqmodel <: DiffEqBase.AbstractEnsembleProblem} <: ScalarMonteCarloProcess{Float64}
	model::absdiffeqmodel
	final_trasform::Function
	underlying::AbstractUnderlying
	function MonteCarloDiffEqModel(model::absdiffeqmodel,final_trasform::Function,underlying::AbstractUnderlying) where {absdiffeqmodel <: DiffEqBase.AbstractEnsembleProblem}
        return new{absdiffeqmodel}(model,final_trasform,underlying)
    end
	function MonteCarloDiffEqModel(model::absdiffeqmodel,underlying::AbstractUnderlying) where {absdiffeqmodel <: DiffEqBase.AbstractEnsembleProblem}
		func(x)=identity(x);
        return new{absdiffeqmodel}(model,func,underlying)
    end
end

export MonteCarloDiffeEqModel;

function simulate!(X,mcProcess::MonteCarloDiffEqModel,rfCurve::ZeroRate,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4,type5},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode, type5 <: Random.AbstractRNG}
	
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	Dt=T/Nstep
	@assert T>0.0
	diffeqmodel=mcProcess.model;
	sol = solve(diffeqmodel,SOSRI(),EnsembleThreads(),trajectories=Nsim,dt=Dt,adaptive=false)
	X.=[path.u[j] for path in sol.u, j in 1:(Nstep+1)];
	f(x)=mcProcess.final_trasform(x);
	broadcast!(f,X,X)
	nothing;
end
