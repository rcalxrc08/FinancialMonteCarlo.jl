
type GeometricBrownianMotion{num,num1<:Number}<:ItoProcess
	sigma::num
	drift::num1
	function GeometricBrownianMotion(sigma::num,drift::num1) where {num,num1 <: Number}
        if sigma <= 0.0
            error("Volatility must be positive")
        else
            return new{num,num1}(sigma,drift)
        end
    end
end

export GeometricBrownianMotion;

function simulate(mcProcess::GeometricBrownianMotion,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::Float64,monteCarloMode::MonteCarloMode=standard)
	if T<=0.0
		error("Final time must be positive");
	end
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	sigma_gbm=mcProcess.sigma;
	mu_gbm=mcProcess.drift;
	drift_bm=mu_gbm-sigma_gbm^2/2;
	const dictGBM=Dict{String,Number}("sigma"=>sigma_gbm, "drift" => drift_bm)
	BrownianData=MonteCarloConfiguration(mcBaseData.Nsim,mcBaseData.Nstep)
	X=simulate(BrownianMotion(sigma_gbm,drift_bm),spotData,BrownianData,T,monteCarloMode)
	S=exp.(X);
	return S;
end