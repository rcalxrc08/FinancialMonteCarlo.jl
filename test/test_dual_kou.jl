using BenchmarkTools, DualNumbers
using MonteCarlo
@show "KouModel"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=dual(0.2,1.0); 
#sigma=0.2; 
p=0.3; 
lam=5.0; 
lamp=30.0; 
lamm=20.0;
const dict=Dict{String,Number}("sigma"=>sigma, "lambda" => lam, "p" => p, "lambdap" => lamp, "lambdam" => lamm)
mc=MonteCarloBaseData(dict,Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);

FwdData=ForwardData(T)
EUData=EUOptionData(T,K)
BarrierData=BarrierOptionData(T,K,D)
AsianData1=AsianFloatingStrikeOptionData(T)
AsianData2=AsianFixedStrikeOptionData(T,K)
Model=KouProcess();

@show FwdPrice=pricer(Model,spotData1,mc,FwdData,Forward());						
@show EuPrice=pricer(Model,spotData1,mc,EUData,EuropeanOption());
@show BarrierPrice=pricer(Model,spotData1,mc,BarrierData,BarrierOption());
@show AsianPrice1=pricer(Model,spotData1,mc,AsianData1,AsianFloatingStrikeOption());
@show AsianPrice2=pricer(Model,spotData1,mc,AsianData2,AsianFixedStrikeOption());

@assert abs(FwdPrice-99.41332633109904)<toll
@assert abs(EuPrice-10.347332240535199)<toll
@assert abs(BarrierPrice-8.860123655599818)<toll
@assert abs(AsianPrice1-5.81798437145069)<toll