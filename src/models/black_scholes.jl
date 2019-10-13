"""
Struct for Black Scholes Process

		bsProcess=BlackScholesProcess(σ::num1) where {num1 <: Number}
	
Where:\n
		σ	=	volatility of the process.
"""
mutable struct BlackScholesProcess{num<:Number}<:ItoProcess
	σ::num
	function BlackScholesProcess(σ::num) where {num <: Number}
        if σ <= 0.0
            error("Volatility must be positive")
        else
            return new{num}(σ)
        end
    end
end

export BlackScholesProcess;

function simulate(mcProcess::BlackScholesProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	if T<=0.0
		error("Final time must be positive");
	end
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ_gbm=mcProcess.σ;
	mu_gbm=r-d;
	
	S=S0.*simulate(GeometricBrownianMotion(σ_gbm,mu_gbm),spotData,mcBaseData,T,parallelMode)
	
	return S;
	
end