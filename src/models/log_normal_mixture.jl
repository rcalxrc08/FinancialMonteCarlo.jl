"""
Struct for LogNormalMixture

		lnmModel=LogNormalMixture(η::Array{num1},λ::Array{num2}) where {num1,num2<: Number}
	
Where:\n
		η  =	Array of volatilities.
		λ  = 	Array of weights.
"""
mutable struct LogNormalMixture{num <: Number,num2 <: Number}<:ItoProcess
	η::Array{num,1}
	λ::Array{num2,1}
	underlying::Underlying
	function LogNormalMixture(η::Array{num,1},λ::Array{num2,1},underlying::Underlying) where {num <: Number,num2 <: Number}
        if minimum(η) <= 0.0
            error("Volatilities must be positive")
        elseif minimum(λ) <= 0.0
            error("weights must be positive")
        elseif sum(λ) > 1.0
            error("λs must be weights")
        elseif length(λ) != length(η)-1
            error("Check vector lengths")
        else
            return new{num,num2}(η,λ,underlying)
        end
    end
	function LogNormalMixture(η::Array{num,1},λ::Array{num2,1},S0::num3) where {num <: Number,num2 <: Number, num3 <:Number}
        if minimum(η) <= 0.0
            error("Volatilities must be positive")
        elseif minimum(λ) <= 0.0
            error("weights must be positive")
        elseif sum(λ) > 1.0
            error("λs must be weights")
        elseif length(λ) != length(η)-1
            error("Check vector lengths")
        else
            return new{num,num2}(η,λ,Underlying(S0))
        end
    end
end

export LogNormalMixture;

function simulate(mcProcess::LogNormalMixture,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	if T<=0.0
		error("Final time must be positive");
	end
	r=spotData.r;
	S0=mcProcess.underlying.S0;
	d=dividend(mcProcess);
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	η_gbm=mcProcess.η;
	λ_gmb=mcProcess.λ
	push!(λ_gmb,1.0-sum(mcProcess.λ))
	mu_gbm=r-d;
	S=λ_gmb[1].*simulate(GeometricBrownianMotion(η_gbm[1],mu_gbm,Underlying(S0)),spotData,mcBaseData,T);
	for (η_gbm_,λ_gmb_) in zip(η_gbm[2:end],λ_gmb[2:end])
		S+=λ_gmb_.*simulate(GeometricBrownianMotion(η_gbm_,mu_gbm,Underlying(S0)),spotData,mcBaseData,T)
	end
	return S;
	
end