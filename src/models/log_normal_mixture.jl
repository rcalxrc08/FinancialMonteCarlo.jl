type LogNormalMixture{num,num2<:Number}<:ItoProcess
	eta::Array{num,1}
	lambda::Array{num2,1}
	function LogNormalMixture(eta::Array{num,1},lambda::Array{num2,1}) where {num,num2 <: Number}
        if minimum(eta) <= 0.0
            error("Volatilities must be positive")
        elseif minimum(lambda) <= 0.0
            error("weights must be positive")
        elseif sum(lambda) > 1.0
            error("lambdas must be weights")
        elseif length(lambda) != length(eta)-1
            error("Check vector lengths")
        else
            return new{num,num2}(eta,lambda)
        end
    end
end

export LogNormalMixture;

function simulate(mcProcess::LogNormalMixture,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::Float64,monteCarloMode::MonteCarloMode=standard)
	if T<=0.0
		error("Final time must be positive");
	end
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	eta_gbm=mcProcess.eta;
	lambda_gmb=mcProcess.lambda
	push!(lambda_gmb,1.0-sum(mcProcess.lambda))
	mu_gbm=r-d;
	S=lambda_gmb[1]*S0.*simulate(GeometricBrownianMotion(eta_gbm[1],mu_gbm),spotData,mcBaseData,T,monteCarloMode);
	for (eta_gbm_,lambda_gmb_) in zip(eta_gbm[2:end],lambda_gmb[2:end])
		S+=lambda_gmb_*S0.*simulate(GeometricBrownianMotion(eta_gbm_,mu_gbm),spotData,mcBaseData,T,monteCarloMode)
	end
	return S;
	
end