using FinancialMonteCarlo,Nabla;
@show "Black Scholes Model"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=10000;
Nstep=30;
sigma=0.2
mc=MonteCarloConfiguration(Nsim,Nstep);
toll=0.8

rfCurve=ZeroRate(S0,r,d);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=BlackScholesProcess(sigma);


#f__(x) = pricer(BlackScholesProcess(x[1]),ZeroRate(x[2],x[3],d),mc,EUData);
f__(x) = pricer(BlackScholesProcess(x[1]),ZeroRate(x[2],x[3],d),mc,FwdData);
#x=Float64[sigma,S0]
x=Float64[sigma,S0,r]
gradient(x -> f__(x), x)
