"""
Struct for Geometric Brownian Motion

		gbmProcess=GeometricBrownianMotion(σ::num1,μ::num2) where {num1,num2 <: Number}
	
Where:\n
		σ	=	volatility of the process.
		μ	=	drift of the process.
"""
mutable struct GeometricBrownianMotion{num <: Number, num1 <: Number, num2 <: Number , num3 <: Number}<:ItoProcess
	σ::num
	μ::num1
	underlying::Underlying{num2,num3}
	function GeometricBrownianMotion(σ::num,μ::num1,underlying::Underlying{num2,num3}) where {num <: Number , num1 <: Number, num2 <: Number, num3 <: Number}
        if σ <= 0.0
            error("Volatility must be positive")
        else
            return new{num,num1,num2,num3}(σ,μ,underlying)
        end
    end
	function GeometricBrownianMotion(σ::num,μ::num1,S0::num2) where {num <: Number , num1 <: Number, num2 <: Number}
        if σ <= 0.0
            error("Volatility must be positive")
        else
            return new{num,num1,num2,Float64}(σ,μ,Underlying(S0))
        end
    end
end

export GeometricBrownianMotion;

function simulate(mcProcess::GeometricBrownianMotion,spotData::ZeroRateCurve,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	if T<=0.0
		error("Final time must be positive");
	end
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ_gbm=mcProcess.σ;
	mu_gbm=mcProcess.μ;
	μ_bm=mu_gbm-σ_gbm^2/2;
	X=simulate(BrownianMotion(σ_gbm,μ_bm,Underlying(0.0)),spotData,mcBaseData,T)
	S=(mcProcess.underlying.S0).*exp.(X);
	return S;
end