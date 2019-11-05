using BenchmarkTools,FinancialMonteCarlo,ReverseDiff

S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=0.2
p=0.3; 
lam=5.0; 
lamp=30.0; 
lamm=20.0;
mc=MonteCarloConfiguration(Nsim,Nstep);
toll=1e-3;

spotData1=equitySpotData(r);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=KouProcess(sigma,lam,p,lamp,lamm,Underlying(S0,d));
Model=KouProcess(x[1],lam,p,lamp,lamm,Underlying(S0,d));

#f(x) = pricer(KouProcess(x[1],lam,p,lamp,lamm,Underlying(S0,d)),equitySpotData(x[2],x[3],x[4]),mc,EuropeanOption(x[5],K));
f(x) = pricer(KouProcess(x[1],lam,p,lamp,lamm,Underlying(S0,d)),equitySpotData(x[2],r,d),mc,EuropeanOption(T,K));
#x=Float64[sigma,S0,r,d,T]
x=Float64[sigma,S0]
g = x -> ReverseDiff.gradient(f, x);
y0=g(x);
@btime f(x);
@btime g(x);

#f1_(x) = pricer(BlackScholesProcess(x[1],Underlying(S0,d)),equitySpotData(x[2],x[3],x[4]),mc,AmericanOption(x[5],K));
f1_(x) = pricer(KouProcess(x[1],lam,p,lamp,lamm,Underlying(S0,d)),equitySpotData(x[2],r,d),mc,AmericanOption(T,K));
g_ = x -> ReverseDiff.gradient(f1_, x);
#y0=g(x);
@btime f1_(x);
#@btime g_(x);