using BenchmarkTools
using MonteCarlo
@show "MertonProcess"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=0.2; 
lam=5.0; 
mu1=0.03; 
sigma1=0.02;
const dict=Dict{String,Number}("sigma"=>sigma, "lambda" => lam, "muJ" => mu1, "sigmaJ" => sigma1)
mc=MonteCarloBaseData(dict,Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);

FwdData=ForwardData(T)
EUData=EUOptionData(T,K)
BarrierData=BarrierOptionData(T,K,D)
AsianData=AsianFloatingStrikeOptionData(T)
Model=MertonProcess();

@show FwdPrice=pricer(Model,spotData1,mc,FwdData,Forward());						
@show EuPrice=pricer(Model,spotData1,mc,EUData,EuropeanOption());
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,BarrierOption());
@show AsianPrice=pricer(Model,spotData1,mc,AsianData,AsianFloatingStrikeOption());

@assert abs(FwdPrice-99.1188767166039)<toll
@assert abs(EuPrice-9.084327245917533)<toll
@assert abs(BarrierPrice-7.880881290426765)<toll
@assert abs(AsianPrice-5.129020349580892)<toll