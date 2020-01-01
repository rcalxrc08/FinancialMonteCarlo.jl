"""
Struct for Geometric Brownian Motion

		gbmProcess=GeometricBrownianMotionVec(σ::num1,μ::num2) where {num1,num2 <: Number}
	
Where:\n
		σ	=	volatility of the process.
		μ	=	drift of the process.
"""
mutable struct GeometricBrownianMotionVec{num <: Number, num1 <: Number, num2 <: Number , num3 <: Number, num4 <: Number}<:ItoProcess
	σ::num
	μ::Curve{num1,num4}
	underlying::Underlying{num2,num3}
	function GeometricBrownianMotionVec(σ::num,μ::Curve{num1,num4},underlying::Underlying{num2,num3}) where {num <: Number , num1 <: Number, num2 <: Number, num3 <: Number, num4 <: Number}
        if σ <= 0.0
            error("Volatility must be positive")
        else
            return new{num,num1,num2,num3,num4}(σ,μ,underlying)
        end
    end
end

function GeometricBrownianMotion(σ::num,μ::Curve{num1,num4},underlying::Underlying{num2,num3}) where {num <: Number, num1 <: Number, num4 <: Number, num2 <: Number , num3 <: Number}
	if σ <= 0.0
		error("Volatility must be positive")
	else
		return GeometricBrownianMotionVec(σ,μ,underlying)
	end
end

export GeometricBrownianMotionVec;

function simulate(mcProcess::GeometricBrownianMotionVec,rfCurve::AbstractZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	if T<=0.0
		error("Final time must be positive");
	end
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ_gbm=mcProcess.σ;
	mu_gbm=mcProcess.μ;
	μ_bm=mu_gbm.-(σ_gbm^2/2);
	X=simulate(BrownianMotion(σ_gbm,μ_bm,Underlying(0.0)),rfCurve,mcBaseData,T)
	S=(mcProcess.underlying.S0).*exp.(X);
	return S;
end