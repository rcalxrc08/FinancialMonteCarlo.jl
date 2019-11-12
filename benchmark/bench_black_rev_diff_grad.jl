using BenchmarkTools,FinancialMonteCarlo,ReverseDiff

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
toll=1e-3;

spotData1=ZeroRateCurve(r);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=BlackScholesProcess(sigma,Underlying(S0,d));

f(x) = pricer(BlackScholesProcess(x[1],Underlying(x[2],d)),ZeroRateCurve(r),mc,EuropeanOption(x[3],K));
#f(x) = pricer(BlackScholesProcess(x[1]),ZeroRateCurve(x[2],r,d),mc,EuropeanOption(T,K),Underlying(S0,d));
#x=Float64[sigma,S0,r,d,T]
x=Float64[sigma,S0,T]
g = x -> ReverseDiff.gradient(f, x);
y0=g(x);
@btime f(x);
@btime g(x);

#f1_(x) = pricer(BlackScholesProcess(x[1]),ZeroRateCurve(x[2],x[3],x[4]),mc,AmericanOption(x[5],K),Underlying(S0,d));
f1_(x) = pricer(BlackScholesProcess(x[1],Underlying(x[2],d)),ZeroRateCurve(r),mc,AmericanOption(x[3],K));
g_ = x -> ReverseDiff.gradient(f1_, x);
#y0=g(x);
@btime f1_(x);
#@btime g_(x);