using BenchmarkTools, DualNumbers
using MonteCarlo

S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=dual(0.2,1.0); 
theta1=0.01; 
k1=0.03; 
sigma1=0.02;
ParamDict=Dict{String,Number}("sigma"=>sigma, "theta" => theta1, "k" => k1)
mc=MonteCarloBaseData(ParamDict,Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);

FwdData=ForwardData(T)
EUData=EUOptionData(T,K)
AMData=AMOptionData(T,K)
BarrierData=BarrierOptionData(T,K,D)
AsianData=AsianFloatingStrikeOptionData(T)
Model=VarianceGammaProcess();

@btime FwdPrice=pricer(Model,spotData1,mc,FwdData,Forward());						
@btime EuPrice=pricer(Model,spotData1,mc,EUData,EuropeanOption());
@btime AmPrice=pricer(Model,spotData1,mc,AMData,AmericanOption());
@btime BarrierPrice=pricer(Model,spotData1,mc,BarrierData,BarrierOption());
@btime AsianPrice=pricer(Model,spotData1,mc,AsianData,AsianFloatingStrikeOption());