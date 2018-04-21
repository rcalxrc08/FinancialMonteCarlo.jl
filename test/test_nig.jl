using BenchmarkTools, MonteCarlo;
@show "NormalInverseGaussianProcess"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=0.2; 
theta1=0.01; 
k1=0.03; 
sigma1=0.02;
ParamDict=Dict{String,Number}("sigma"=>sigma, "theta" => theta1, "k" => k1)
mc=MonteCarloBaseData(ParamDict,Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);

FwdData=ForwardData(T)
EUData=EUOptionData(T,K)
BarrierData=BarrierOptionData(T,K,D)
AsianData=AsianFloatingStrikeOptionData(T)
Model=NormalInverseGaussianProcess();

@show FwdPrice=pricer(Model,spotData1,mc,FwdData,Forward());						
@show EuPrice=pricer(Model,spotData1,mc,EUData,EuropeanOption());
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,BarrierOption());
@show AsianPrice=pricer(Model,spotData1,mc,AsianData,AsianFloatingStrikeOption());

@assert abs(FwdPrice-97.95852227697686)<toll
@assert abs(EuPrice-7.738298817933206)<toll
@assert abs(BarrierPrice-6.886023820038332)<toll
@assert abs(AsianPrice-4.414948846776423)<toll