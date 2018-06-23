function runnerMonteCarlo(Model::AbstractMonteCarloProcess,mc::MonteCarloBaseData)
	#@show Model
	S0=100.0;
	K=100.0;
	r=0.02;
	T=1.0;
	
	d=0.01;
	
	D=90.0;
	spotData1=equitySpotData(S0,r,d);

	FwdData=Forward(T)
	EUData=EuropeanOption(T,K)
	AMData=AmericanOption(T,K)
	BarrierData=BarrierOptionDownOut(T,K,D)
	AsianData=AsianFloatingStrikeOption(T)

	optionDatas=[FwdData,EUData,AMData,BarrierData,AsianData]
	
	(FwdPrice,EuPrice,AMPrice,BarrierPrice,AsianPrice)=pricer(Model,spotData1,mc,optionDatas)

	@show FwdPrice
	@show EuPrice
	@show AMPrice
	@show BarrierPrice
	@show AsianPrice

end

using DualNumbers;

function runnerMonteCarloDual(Model::AbstractMonteCarloProcess,mc::MonteCarloBaseData)
	#@show Model
	S0=dual(100.0,1.0);
	K=100.0;
	r=0.02;
	T=1.0;
	
	d=0.01;
	
	D=90.0;
	spotData1=equitySpotData(S0,r,d);

	FwdData=Forward(T)
	EUData=EuropeanOption(T,K)
	AMData=AmericanOption(T,K)
	BarrierData=BarrierOptionDownOut(T,K,D)
	AsianData=AsianFloatingStrikeOption(T)

	optionDatas=[FwdData,EUData,AMData,BarrierData,AsianData]
	
	(FwdPrice,EuPrice,AMPrice,BarrierPrice,AsianPrice)=pricer(Model,spotData1,mc,optionDatas)
	

	@show FwdPrice
	@show EuPrice
	@show AMPrice
	@show BarrierPrice
	@show AsianPrice

end