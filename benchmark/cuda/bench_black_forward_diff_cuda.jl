using BenchmarkTools,FinancialMonteCarlo,CuArrays,ForwardDiff
CuArrays.allowscalar(false)
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=100000;
Nstep=30;
sigma=0.2
mc=MonteCarloConfiguration(Nsim,Nstep);
mc_1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.CudaMode());
mc_2=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.CudaMode_2());
toll=1e-3;

spotData1=equitySpotData(r);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=BlackScholesProcess(sigma,Underlying(S0,d));

f(x::Vector) = pricer(BlackScholesProcess(x[1]),equitySpotData(x[2],x[3],x[4]),mc,EuropeanOption(x[5],K),FinancialMonteCarlo.CudaMode_2());
x=Float64[sigma,S0,r,d,T]
g = x -> ForwardDiff.gradient(f, x);
y0=g(x);
@btime f(x);
@btime g(x);

f_(x::Vector) = pricer(BlackScholesProcess(x[1]),equitySpotData(x[2],x[3],x[4]),mc,AmericanOption(x[5],K),FinancialMonteCarlo.CudaMode_2());
g_ = x -> ForwardDiff.gradient(f_, x);
y0=g_(x);
@btime f_(x);
@btime g_(x);
