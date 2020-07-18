using BenchmarkTools,FinancialMonteCarlo,DualNumbers

S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=dual(0.2,1.0);
mc=MonteCarloConfiguration(Nsim,Nstep);
toll=1e-3;

rfCurve=ZeroRate(r);

FwdData=Forward(T*2.0)
EUData=EuropeanOption(T*1.1,K)
AMData=AmericanOption(T*4.0,K)
BarrierData=BarrierOptionDownOut(T*0.5,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T*3.3)
AsianFixedStrikeData=AsianFixedStrikeOption(T*1.2,K)
Model=BlackScholesProcess(sigma,Underlying(S0,d));

optionDatas=[FwdData,EUData,AMData,BarrierData,AsianFloatingStrikeData,AsianFixedStrikeData]

@btime (FwdPrice,EuPrice,AMPrice,BarrierPrice,AsianPrice1,AsianPrice2)=pricer(Model,rfCurve,mc,optionDatas);