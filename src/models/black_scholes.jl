"""
Struct for Black Scholes Process

		bsProcess=BlackScholesProcess(σ::num1) where {num1 <: Number}
	
Where:\n
		σ	=	volatility of the process.
"""
mutable struct BlackScholesProcess{num <: Number, abstrUnderlying <: AbstractUnderlying} <: ItoProcess{num}
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

function simulate!(S,mcProcess::BlackScholesProcess,rfCurve::AbstractZeroRateCurve,mcBaseData::AbstractMonteCarloConfiguration,T::Number)
	@assert T>0.0
	r=rfCurve.r;
	d=dividend(mcProcess);
	σ_gbm=mcProcess.σ;
	mu_gbm=r-d;
	
	simulate!(S,GeometricBrownianMotion(σ_gbm,mu_gbm,mcProcess.underlying.S0),mcBaseData,T)
	
	nothing
	
end