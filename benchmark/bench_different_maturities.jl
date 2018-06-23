using BenchmarkTools,MonteCarlo,DualNumbers

S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=dual(0.2,1.0);
ParamDict=Dict{String,Number}("sigma"=>sigma)
mc=MonteCarloBaseData(ParamDict,Nsim,Nstep);
toll=1e-3;

spotData1=equitySpotData(S0,r,d);

FwdData=ForwardData(T*2.0)
EUData=EUOptionData(T*1.1,K)
AMData=AmericanStdOption(T*4.0,K)
BarrierData=BarrierOptionData(T*0.5,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOptionData(T*3.3)
AsianFixedStrikeData=AsianFixedStrikeOptionData(T*1.2,K)
Model=BlackScholesProcess();

@show FwdPrice=pricer(Model,spotData1,mc,FwdData,Forward());						
@show EuPrice=pricer(Model,spotData1,mc,EUData);
@show AmPrice=pricer(Model,spotData1,mc,AMData);
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData);
@show AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeData);
@show AsianPrice2=pricer(Model,spotData1,mc,AsianFixedStrikeData);

optionDatas=[FwdData,EUData,AMData,BarrierData,AsianFloatingStrikeData,AsianFixedStrikeData]

(FwdPrice,EuPrice,AMPrice,BarrierPrice,AsianPrice1,AsianPrice2)=pricer(Model,spotData1,mc,optionDatas)
@btime (FwdPrice,EuPrice,AMPrice,BarrierPrice,AsianPrice1,AsianPrice2)=pricer(Model,spotData1,mc,optionDatas)
@show FwdPrice
@show EuPrice
@show AmPrice
@show BarrierPrice
@show AsianPrice1
@show AsianPrice2