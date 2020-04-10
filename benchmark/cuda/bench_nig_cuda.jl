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
theta1=0.01; 
k1=0.03; 
sigma1=0.02;
mc=MonteCarloConfiguration(Nsim,Nstep);
toll=0.8;
mc=MonteCarloConfiguration(Nsim,Nstep);

mc_2=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.CudaMode_2());
rfCurve=ZeroRate(r);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=NormalInverseGaussianProcess(sigma,theta1,k1,Underlying(S0,d));

@show "CUDA_1 fwd"

@show "STD fwd"
@btime FwdPrice=pricer(Model,rfCurve,mc,FwdData);
@show "CUDA_1 fwd"

@show "CUDA_2 fwd"
@btime FwdPrice=pricer(Model,rfCurve,mc_2,FwdData);
@show "std eu"
@btime EuPrice=pricer(Model,rfCurve,mc,EUData);
@show "CUDA_1 eu"

@show "CUDA_2 eu"
@btime EuPrice=pricer(Model,rfCurve,mc_2,EUData);
@show "std am"
@btime AmPrice=pricer(Model,rfCurve,mc,AMData);
@show "CUDA_1 am"

@show "CUDA_2 am"
@btime AmPrice=pricer(Model,rfCurve,mc_2,AMData);