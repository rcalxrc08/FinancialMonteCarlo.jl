using BenchmarkTools
using MonteCarlo

S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=100000;
Nstep=30;
sigma=0.2;
mc=MonteCarloConfiguration(Nsim,Nstep);
toll=1e-3;

spotData1=equitySpotData(S0,r,d);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=BlackScholesProcess(sigma);

@btime FwdPrice=pricer(Model,spotData1,mc,FwdData,MonteCarlo.parallel_cuda_gpu);

@btime FwdPrice=pricer(Model,spotData1,mc,FwdData);
@btime FwdPrice=pricer(Model,spotData1,mc,FwdData,MonteCarlo.parallel_cuda_gpu);
@btime FwdPrice=pricer(Model,spotData1,mc,FwdData,MonteCarlo.parallel_cuda_gpu_trial);
@btime EuPrice=pricer(Model,spotData1,mc,EUData);
@btime EuPrice=pricer(Model,spotData1,mc,EUData,MonteCarlo.parallel_cuda_gpu);
@btime EuPrice=pricer(Model,spotData1,mc,EUData,MonteCarlo.parallel_cuda_gpu_trial);
@btime AmPrice=pricer(Model,spotData1,mc,AMData);
@btime AmPrice=pricer(Model,spotData1,mc,AMData,MonteCarlo.parallel_cuda_gpu);
@btime AmPrice=pricer(Model,spotData1,mc,AMData,MonteCarlo.parallel_cuda_gpu_trial);