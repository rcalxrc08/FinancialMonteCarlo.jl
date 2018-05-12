type AmericanOption<:PayoffMC end

struct AMOptionData<:AbstractEuropeanOptionData
	T::Float64
	K::Float64
	r::Float64
end

export AmericanOption,AMOptionData;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,amOptionData,AmericanOption,isCall=true)
	
Where:\n
		S           = Paths of the Underlying.
		amOptionData  = Datas of the Option.
		AmericanOption = Type of the Option
		isCall = true for Call Options, false for price Options.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},amOptionData::AMOptionData,Payoff::AmericanOption,isCall::Bool=true) where{num<:Number}
	iscall=isCall?1:-1
	S0=S[1,1];
	Nsim=length(S[1:end,1]);
	Nstep=length(S[1,1:end])
	K=amOptionData.K;
	r=amOptionData.r;
	dt=amOptionData.T/Nstep
	# INIZIALIZZAZIONE ========================
	Istanti_esercizio=Nstep*ones(Nsim,1);
	price=max.(0.0,iscall*(S[:,end]-K)); #payoff
	

	# PROCEDURA BACKWARD ======================
	for j in Nstep-1:-1:1
		Inmoney=find((S[:,j]-K)*iscall.>0.0);
		if Inmoney!=[]
			S_I=S[Inmoney,j];
			#-- Intrinsic Value
			IV=(S_I-K)*iscall;
			#-- Continuation Value 
			#- Regressione Lineare
			A=[ones(length(S_I),1) S_I S_I.^2];
			b=price[Inmoney].*exp.(-r*dt*(Istanti_esercizio[Inmoney]-j));
			MAT=A'*A;			
			LuMat=lufact(MAT);
			Btilde=A'*b;
			alpha=LuMat\Btilde;
			#alpha=A\b;
			#- Valore di continuazione
			CV=A*alpha;
			#----------
			# istante j è di esercizio anticipato?
			Index=find(IV.>CV);
			Esercizio_Anticipato=Inmoney[Index];
			# Update
			price[Esercizio_Anticipato]=IV[Index];
			Istanti_esercizio[Esercizio_Anticipato]=j;
		end
	end
	prezzo=max.(iscall*(S0-K),mean(price.*exp.(-r*dt*Istanti_esercizio)))


	
	return prezzo;
end