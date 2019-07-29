using BenchmarkTools
using FinancialMonteCarlo,CuArrays
CuArrays.allowscalar(false)

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
mc=MonteCarloConfiguration(Nsim,Nstep);
toll=0.8;

spotData1=equitySpotData(S0,r,d);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=MertonProcess(sigma,lam,mu1,sigma1);


@btime FwdPrice=pricer(Model,spotData1,mc,FwdData);
@btime FwdPrice=pricer(Model,spotData1,mc,FwdData,FinancialMonteCarlo.CudaMode());
@btime FwdPrice=pricer(Model,spotData1,mc,FwdData,FinancialMonteCarlo.CudaMode_2());
@btime EuPrice=pricer(Model,spotData1,mc,EUData);
@btime EuPrice=pricer(Model,spotData1,mc,EUData,FinancialMonteCarlo.CudaMode());
@btime EuPrice=pricer(Model,spotData1,mc,EUData,FinancialMonteCarlo.CudaMode_2());
@btime AmPrice=pricer(Model,spotData1,mc,AMData);
@btime AmPrice=pricer(Model,spotData1,mc,AMData,FinancialMonteCarlo.CudaMode());