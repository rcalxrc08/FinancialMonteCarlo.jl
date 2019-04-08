using BenchmarkTools
using FinancialMonteCarlo,CuArrays
CuArrays.allowscalar(false)
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=1000000;
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
@show "CUDA_1 fwd"
@btime FwdPrice=pricer(Model,spotData1,mc,FwdData,FinancialMonteCarlo.standard,FinancialMonteCarlo.CudaMode());

@show "STD fwd"
@btime FwdPrice=pricer(Model,spotData1,mc,FwdData);
@show "CUDA_1 fwd"
@btime FwdPrice=pricer(Model,spotData1,mc,FwdData,FinancialMonteCarlo.standard,FinancialMonteCarlo.CudaMode());
@show "CUDA_2 fwd"
@btime FwdPrice=pricer(Model,spotData1,mc,FwdData,FinancialMonteCarlo.standard,FinancialMonteCarlo.CudaMode_2());
@show "std eu"
@btime EuPrice=pricer(Model,spotData1,mc,EUData);
@show "CUDA_1 eu"
@btime EuPrice=pricer(Model,spotData1,mc,EUData,FinancialMonteCarlo.standard,FinancialMonteCarlo.CudaMode());
@show "CUDA_2 eu"
@btime EuPrice=pricer(Model,spotData1,mc,EUData,FinancialMonteCarlo.standard,FinancialMonteCarlo.CudaMode_2());
@show "std am"
@btime AmPrice=pricer(Model,spotData1,mc,AMData);
@show "CUDA_1 am"
@btime AmPrice=pricer(Model,spotData1,mc,AMData,FinancialMonteCarlo.standard,FinancialMonteCarlo.CudaMode());
@show "CUDA_2 am"
@btime AmPrice=pricer(Model,spotData1,mc,AMData,FinancialMonteCarlo.standard,FinancialMonteCarlo.CudaMode_2());