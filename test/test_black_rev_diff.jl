using FinancialMonteCarlo,ReverseDiff;
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

spotData1=equitySpotData(S0,r,d);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
AMData=AmericanOption(T,K)
BarrierData=BarrierOptionDownOut(T,K,D)
AsianFloatingStrikeData=AsianFloatingStrikeOption(T)
AsianFixedStrikeData=AsianFixedStrikeOption(T,K)
Model=BlackScholesProcess(sigma);

x=Float64[sigma,S0,r,d]

@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1]),equitySpotData(x[2],x[3],x[4]),mc,FwdData), x)
@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1]),equitySpotData(x[2],x[3],x[4]),mc,EUData), x)
@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1]),equitySpotData(x[2],x[3],x[4]),mc,AMData), x)
@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1]),equitySpotData(x[2],x[3],x[4]),mc,BarrierData), x)
@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1]),equitySpotData(x[2],x[3],x[4]),mc,AsianFloatingStrikeData), x)
@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1]),equitySpotData(x[2],x[3],x[4]),mc,AsianFixedStrikeData), x)