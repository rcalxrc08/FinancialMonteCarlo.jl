"""
Struct for Asian Fixed Strike Option

		asOption=AsianFixedStrikeOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
mutable struct AsianFixedStrikeOption{T1,T2<:Number}<:AsianPayoff
	T::T1
	K::T2
	isCall::Bool
	function AsianFixedStrikeOption(T::T1,K::T2,isCall::Bool=true) where {T1, T2 <: Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
		
            return new{T1,T2}(T,K,isCall)
        end
    end

end

export AsianFixedStrikeOption;


function payoff(S::AbstractMatrix{num},asianFixedStrikePayoff::AsianFixedStrikeOption,spotData::ZeroRateCurve,T1::num2=asianFixedStrikePayoff.T) where{num <: Number,num2 <: Number}
	iscall=asianFixedStrikePayoff.isCall ? 1 : -1
	r=spotData.r;
	T=asianFixedStrikePayoff.T;
	K=asianFixedStrikePayoff.K;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	@inbounds f(S::abstractArray) where {abstractArray<:AbstractArray{num_}} where {num_<:Number}=max(iscall*(mean(S)-K),0.0)
	@inbounds payoff2=[f(view(S,i,1:index1)) for i in 1:Nsim];
	
	return payoff2*exp(-r*T);
end