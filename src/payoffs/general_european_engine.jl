
function payoff(S::AbstractMatrix{num},euPayoff::EuropeanPayoff,spotData::ZeroRateCurve,T1::num2=euPayoff.T) where{ num <: Number, num2 <: Number}
	r=spotData.r;
	T=euPayoff.T;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(UInt,T/T1 * NStep)+1;
	@views ST=S[:,index1];
	phi(Sti::numtype_) where {numtype_<:Number}=payout(Sti,euPayoff);
	payoff2=phi.(ST);
	
	return payoff2*exp(-r*T);
end