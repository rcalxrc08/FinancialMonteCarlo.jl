using BenchmarkTools,MonteCarlo,DualNumbers;

S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=0.2; 
sigma_zero=dual(0.2,1.0);
kappa=0.01;
theta=0.03;
lambda=0.01;
rho=0.0;
ParamDict=Dict{String,Number}("sigma"=>sigma,"theta"=>theta, "kappa"=>kappa, "sigma_zero" => sigma_zero, "lambda" => lambda, "rho" => rho)
mc=MonteCarloBaseData(ParamDict,Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);

FwdData=ForwardData(T)
EUData=EUOptionData(T,K)
AMData=AMOptionData(T,K)
BarrierData=BarrierOptionData(T,K,D)
AsianData1=AsianFloatingStrikeOptionData(T)
AsianData2=AsianFixedStrikeOptionData(T,K)
Model=HestonProcess();

@btime FwdPrice=pricer(Model,spotData1,mc,FwdData,Forward());						
@btime EuPrice=pricer(Model,spotData1,mc,EUData,EuropeanOption());
@btime AmPrice=pricer(Model,spotData1,mc,AMData,AmericanOption());
@btime BarrierPrice=pricer(Model,spotData1,mc,BarrierData,BarrierOptionDownOut());
@btime AsianPrice1=pricer(Model,spotData1,mc,AsianData1,AsianFloatingStrikeOption());