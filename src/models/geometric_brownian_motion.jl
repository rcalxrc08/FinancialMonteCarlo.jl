"""
Struct for Geometric Brownian Motion

		gbmProcess=GeometricBrownianMotion(σ::num1,μ::num2) where {num1,num2 <: Number}
	
Where:\n
		σ	=	volatility of the process.
		μ	=	drift of the process.
"""
mutable struct GeometricBrownianMotion{num <: Number, num1 <: Number, abstrUnderlying <: AbstractUnderlying}<:ItoProcess
	σ::num
	μ::num1
	underlying::abstrUnderlying
	function GeometricBrownianMotion(σ::num,μ::num1,underlying::abstrUnderlying) where {num <: Number , num1 <: Number, abstrUnderlying <: AbstractUnderlying}
        if σ <= 0.0
            error("Volatility must be positive")
        else
            return new{num,num1,abstrUnderlying}(σ,μ,underlying)
        end
    end
end

export GeometricBrownianMotion;

function simulate(mcProcess::GeometricBrownianMotion,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4,type5},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode, type5 <: Random.AbstractRNG}
	if T<=0.0
		error("Final time must be positive");
	end
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ_gbm=mcProcess.σ;
	mu_gbm=mcProcess.μ;
	μ_bm=mu_gbm-σ_gbm^2/2;
	X=simulate(BrownianMotion(σ_gbm,μ_bm,Underlying(0.0)),mcBaseData,T)
	S=(mcProcess.underlying.S0).*exp.(X);
	return S;
end