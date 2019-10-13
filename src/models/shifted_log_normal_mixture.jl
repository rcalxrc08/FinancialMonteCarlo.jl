"""
Struct for ShiftedLogNormalMixture

		lnmModel=ShiftedLogNormalMixture(η::Array{num1},λ::Array{num2},num3) where {num1,num2,num3<: Number}
	
Where:\n
		η  =	Array of volatilities.
		λ  = 	Array of weights.
		α  = 	shift.
"""
mutable struct ShiftedLogNormalMixture{num <: Number,num2 <: Number, num3 <: Number}<:ItoProcess
	η::Array{num,1}
	λ::Array{num2,1}
	α::num3
	function ShiftedLogNormalMixture(η::Array{num,1},λ::Array{num2,1},α::num3) where {num <: Number,num2 <: Number, num3 <: Number}
        if minimum(η) <= 0.0
            error("Volatilities must be positive")
        elseif minimum(λ) <= 0.0
            error("weights must be positive")
        elseif sum(λ) > 1.0
            error("λs must be weights")
        elseif length(λ) != length(η)-1
            error("Check vector lengths")
        else
            return new{num,num2,num3}(η,λ,α)
        end
    end
end

export ShiftedLogNormalMixture;

function simulate(mcProcess::ShiftedLogNormalMixture,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode}
	if T<=0.0
		error("Final time must be positive");
	end
	r=spotData.r;
	S0=spotData.S0;
	d=spotData.d;
	Nsim=mcBaseData.Nsim;
	Nstep=mcBaseData.Nstep;
	η_gbm=copy(mcProcess.η);
	λ_gmb=copy(mcProcess.λ)
	push!(λ_gmb,1.0-sum(mcProcess.λ))
	mu_gbm=r-d;
	dt=T/Nstep
	A0=S0*(1-mcProcess.α)
	S=(λ_gmb[1]*A0).*simulate(GeometricBrownianMotion(η_gbm[1],mu_gbm),equitySpotData(A0,r,d),mcBaseData,T);
	for (η_gbm_,λ_gmb_) in zip(η_gbm[2:end],λ_gmb[2:end])
		S+=(λ_gmb_*A0).*simulate(GeometricBrownianMotion(η_gbm_,mu_gbm),equitySpotData(A0,r,d),mcBaseData,T)
	end
	return S.+mcProcess.α.*S0.*exp.(mu_gbm.*(0.0:dt:T))';
	
end