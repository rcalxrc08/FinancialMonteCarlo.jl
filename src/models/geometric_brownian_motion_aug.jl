"""
Struct for Geometric Brownian Motion

		gbmProcess=GeometricBrownianMotionVec(σ::num1,μ::num2) where {num1,num2 <: Number}
	
Where:\n
		σ	=	volatility of the process.
		μ	=	drift of the process.
"""
mutable struct GeometricBrownianMotionVec{num <: Number, num1 <: Number, num4 <: Number,  abstrUnderlying <: AbstractUnderlying}<:ItoProcess
	σ::num
	μ::Curve{num1,num4}
	underlying::abstrUnderlying
	function GeometricBrownianMotionVec(σ::num,μ::Curve{num1,num4},underlying::abstrUnderlying) where {num <: Number , num1 <: Number, num4 <: Number, abstrUnderlying <: AbstractUnderlying}
        if σ <= 0.0
            error("Volatility must be positive")
        else
            return new{num,num1,num4,abstrUnderlying}(σ,μ,underlying)
        end
    end
end

function GeometricBrownianMotion(σ::num,μ::Curve{num1,num4},underlying::abstrUnderlying) where {num <: Number, num1 <: Number, num4 <: Number, abstrUnderlying <: AbstractUnderlying}
	return GeometricBrownianMotionVec(σ,μ,underlying)
end

export GeometricBrownianMotionVec;

function simulate!(X,mcProcess::GeometricBrownianMotionVec,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4,type5},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode, type5 <: Random.AbstractRNG}
	@assert T>0.0
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ_gbm=mcProcess.σ;
	mu_gbm=mcProcess.μ;
	μ_bm=mu_gbm-(σ_gbm^2/2);
	simulate!(X,BrownianMotion(σ_gbm,μ_bm,Underlying(0.0)),mcBaseData,T)
	S0=mcProcess.underlying.S0;
	f(x)=S0*exp(x);
	broadcast!(f,X,X)
	nothing;
end