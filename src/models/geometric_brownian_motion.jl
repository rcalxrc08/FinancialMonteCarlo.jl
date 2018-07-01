
type GeometricBrownianMotion{num,num1<:Number}<:ItoProcess
	σ::num
	μ::num1
	function GeometricBrownianMotion(σ::num,μ::num1) where {num,num1 <: Number}
        if σ <= 0.0
            error("Volatility must be positive")
        else
            return new{num,num1}(σ,μ)
        end
    end
end

export GeometricBrownianMotion;

function simulate(mcProcess::GeometricBrownianMotion,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::numb,monteCarloMode::MonteCarloMode=standard) where {numb<:Number}
	if T<=0.0
		error("Final time must be positive");
	end
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	σ_gbm=mcProcess.σ;
	mu_gbm=mcProcess.μ;
	μ_bm=mu_gbm-σ_gbm^2/2;
	X=simulate(BrownianMotion(σ_gbm,μ_bm),spotData,mcBaseData,T,monteCarloMode)
	S=exp.(X);
	return S;
end