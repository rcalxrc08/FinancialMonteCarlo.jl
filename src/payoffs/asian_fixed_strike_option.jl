"""
Class for Dispatching Forward Payoff

		forward=Forward(T::num) where {num<:Number}
	
Where:\n
		T	=	Time to maturity of the Forward.
"""
struct AsianFixedStrikeOption{T1,T2<:Number}<:AsianPayoff
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

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,asianFixedStrikePayoff,spotData,T1)
	
Where:\n
		S           = Paths of the Underlying.
		asianFixedStrikePayoff  = Datas of the Option.
		spotData  = Datas of the Spot.
		T1=Final Time of Spot Simulation (default equals Time to Maturity of the option)

		Payoff      = payoff of the Option.
```
"""
function payoff(S::AbstractMatrix{num},asianFixedStrikePayoff::AsianFixedStrikeOption,spotData::equitySpotData,T1::num2=asianFixedStrikePayoff.T) where{num,num2<:Number}
	iscall=asianFixedStrikePayoff.isCall ? 1 : -1
	r=spotData.r;
	T=asianFixedStrikePayoff.T;
	K=asianFixedStrikePayoff.K;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	@inbounds f(S::Array{num})::num=(iscall*(mean(S)-K)>0.0) ? (iscall*(mean(S)-K)) : 0.0;		
	@inbounds payoff2=[f(S1[i,:]) for i in 1:Nsim];
	
	return payoff2*exp(-r*T);
end