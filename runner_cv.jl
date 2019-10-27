using BenchmarkTools
using FinancialMonteCarlo

S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=100000;
Nstep=30;
sigma=0.2;
variate_=FinancialMonteCarlo.ControlVariates(Forward(T),MonteCarloConfiguration(1000,100))
mc=MonteCarloConfiguration(Nsim,Nstep,variate_);
toll=1e-3;

spotData1=equitySpotData(S0,r,d);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=BlackScholesProcess(sigma);

 EuPrice=pricer(Model,spotData1,mc,EUData);
@btime AmPrice=pricer(Model,spotData1,mc,AMData);
@btime BarrierPrice=pricer(Model,spotData1,mc,BarrierData);
@btime AsianPrice1=pricer(Model,spotData1,mc,AsianFloatingStrikeData);
@btime AsianPrice1=pricer(Model,spotData1,mc,AsianFixedStrikeData);

optionDatas=[FwdData,EUData,AMData,BarrierData,AsianFloatingStrikeData,AsianFixedStrikeData]

@btime (FwdPrice,EuPrice,AMPrice,BarrierPrice,AsianPrice1,AsianPrice2)=pricer(Model,spotData1,mc,optionDatas)