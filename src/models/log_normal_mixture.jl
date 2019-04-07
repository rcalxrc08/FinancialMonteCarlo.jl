mutable struct LogNormalMixture{num,num2<:Number}<:ItoProcess
	η::Array{num,1}
	λ::Array{num2,1}
	function LogNormalMixture(η::Array{num,1},λ::Array{num2,1}) where {num,num2 <: Number}
        if minimum(η) <= 0.0
            error("Volatilities must be positive")
        elseif minimum(λ) <= 0.0
            error("weights must be positive")
        elseif sum(λ) > 1.0
            error("λs must be weights")
        elseif length(λ) != length(η)-1
            error("Check vector lengths")
        else
            return new{num,num2}(η,λ)
        end
    end
end

export LogNormalMixture;

function simulate(mcProcess::LogNormalMixture,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::numb,monteCarloMode::MonteCarloMode=standard,parallelMode::BaseMode=SerialMode()) where {numb<:Number}
	if T<=0.0
		error("Final time must be positive");
	end
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	η_gbm=mcProcess.η;
	λ_gmb=mcProcess.λ
	push!(λ_gmb,1.0-sum(mcProcess.λ))
	mu_gbm=r-d;
	S=λ_gmb[1]*S0.*simulate(GeometricBrownianMotion(η_gbm[1],mu_gbm),spotData,mcBaseData,T,monteCarloMode);
	for (η_gbm_,λ_gmb_) in zip(η_gbm[2:end],λ_gmb[2:end])
		S+=λ_gmb_*S0.*simulate(GeometricBrownianMotion(η_gbm_,mu_gbm),spotData,mcBaseData,T,monteCarloMode)
	end
	return S;
	
end