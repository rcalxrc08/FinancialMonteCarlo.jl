"""
Struct for Black Scholes Process

		bsProcess=BlackScholesProcess(σ::num1) where {num1 <: Number}
	
Where:\n
		σ	=	volatility of the process.
"""
mutable struct BlackScholesProcess{num <: Number, abstrUnderlying <: AbstractUnderlying}<:ItoProcess
	σ::num
	underlying::abstrUnderlying
	function BlackScholesProcess(σ::num,underlying::abstrUnderlying) where {num <: Number,  abstrUnderlying <: AbstractUnderlying}
        if σ <= 0.0
            error("Volatility must be positive")
        else
            return new{num,abstrUnderlying}(σ,underlying)
        end
    end
end

export BlackScholesProcess;

function simulate!(S,mcProcess::BlackScholesProcess,rfCurve::AbstractZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4,type5},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode, type5 <: Random.AbstractRNG}
	@assert T>0.0
	r=rfCurve.r;
	d=dividend(mcProcess);
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ_gbm=mcProcess.σ;
	mu_gbm=r-d;
	
	simulate!(S,GeometricBrownianMotion(σ_gbm,mu_gbm,mcProcess.underlying),mcBaseData,T)
	
	nothing
	
end