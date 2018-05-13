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
p=0.3; 
lam=5.0; 
lamp=30.0; 
lamm=20.0;
ParamDict=Dict{String,Number}("sigma"=>sigma, "lambda" => lam, "p" => p, "lambdap" => lamp, "lambdam" => lamm)
mc=MonteCarloBaseData(ParamDict,Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);

FwdData=ForwardData(T)
EUData=EUOptionData(T,K)
AMData=AMOptionData(T,K)
BarrierData=BarrierOptionData(T,K,D)
AsianData1=AsianFloatingStrikeOptionData(T)
AsianData2=AsianFixedStrikeOptionData(T,K)
Model=KouProcess();

@btime FwdPrice=pricer(Model,spotData1,mc,FwdData,Forward());						
@btime EuPrice=pricer(Model,spotData1,mc,EUData,EuropeanOption());
@btime AmPrice=pricer(Model,spotData1,mc,AMData,AmericanOption());
@btime BarrierPrice=pricer(Model,spotData1,mc,BarrierData,BarrierOption());
@btime AsianPrice1=pricer(Model,spotData1,mc,AsianData1,AsianFloatingStrikeOption());