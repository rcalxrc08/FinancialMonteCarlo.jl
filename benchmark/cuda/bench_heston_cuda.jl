using BenchmarkTools,FinancialMonteCarlo,DualNumbers,CuArrays;
CuArrays.allowscalar(false)
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=100000;
Nstep=30;
sigma=0.2; 
sigma_zero=0.2;
kappa=0.01;
theta=0.03;
lambda=0.01;
rho=0.0;
mc=MonteCarloConfiguration(Nsim,Nstep);
mc_1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.CudaMode());
mc_2=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.CudaMode_2());
toll=0.8;

spotData1=equitySpotData(r);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=HestonProcess(sigma,sigma_zero,lambda,kappa,rho,theta,Underlying(S0,d));

@show "CUDA_1 fwd"
@btime FwdPrice=pricer(Model,spotData1,mc_1,FwdData);
@show "STD fwd"
@btime FwdPrice=pricer(Model,spotData1,mc,FwdData);
@show "CUDA_1 fwd"
@btime FwdPrice=pricer(Model,spotData1,mc_1,FwdData);
@show "CUDA_2 fwd"
@btime FwdPrice=pricer(Model,spotData1,mc_2,FwdData);
@show "std eu"
@btime EuPrice=pricer(Model,spotData1,mc,EUData);
@show "CUDA_1 eu"
@btime EuPrice=pricer(Model,spotData1,mc_1,EUData);
@show "CUDA_2 eu"
@btime EuPrice=pricer(Model,spotData1,mc_2,EUData);
@show "std am"
@btime AmPrice=pricer(Model,spotData1,mc,AMData);
@show "CUDA_1 am"
@btime AmPrice=pricer(Model,spotData1,mc_1,AMData);
@show "CUDA_2 am"
@btime AmPrice=pricer(Model,spotData1,mc_2,AMData);