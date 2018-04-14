using BenchmarkTools
using MonteCarlo

S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=0.2;
ParamDict=Dict{String,Number}("sigma"=>sigma)
mc=MonteCarloBaseData(dict,Nsim,Nstep);
toll=1e-3;

spotData1=equitySpotData(S0,r,d);

FwdData=ForwardData(T)
EUData=EUOptionData(T,K)
BarrierData=BarrierOptionData(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOptionData(T)
AsianFixedStrikeData=AsianFixedStrikeOptionData(T,K)
Model=BlackScholesProcess();

@btime FwdPrice=pricer(Model,spotData1,mc,FwdData,Forward());						
@btime EuPrice=pricer(Model,spotData1,mc,EUData,EuropeanOption());
@btime BarrierPrice=pricer(Model,spotData1,mc,BarrierData,BarrierOption());
@btime AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeData,AsianFloatingStrikeOption());
@btime AsianPrice2=pricer(Model,spotData1,mc,AsianFloatingStrikeData,AsianFloatingStrikeOption());