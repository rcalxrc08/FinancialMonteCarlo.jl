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

function simulate!(X,mcProcess::GeometricBrownianMotion,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4,type5},T::numb) where {numb <: Number, type1 <: Integer, type2<: Integer, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode, type5 <: Random.AbstractRNG}
	@assert T>0.0
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ_gbm=mcProcess.σ;
	mu_gbm=mcProcess.μ;
	μ_bm=mu_gbm-σ_gbm^2/2;
	simulate!(X,BrownianMotion(σ_gbm,μ_bm,Underlying(0.0)),mcBaseData,T)
	S0=mcProcess.underlying.S0;
	f(x)=S0*exp(x);
	broadcast!(f,X,X)
	
	nothing;
end