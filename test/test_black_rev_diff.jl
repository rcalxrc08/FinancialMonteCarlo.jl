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
mc1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.AntitheticMC());
toll=0.8

spotData1=equitySpotData(r,d);

Model=BlackScholesProcess(sigma,S0);

x=Float64[sigma,S0,r,d,T]

@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1]),equitySpotData(x[2],x[3],x[4]),mc,Forward(x[5])), x)
@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1]),equitySpotData(x[2],x[3],x[4]),mc,EuropeanOption(x[5],K)), x)
@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1]),equitySpotData(x[2],x[3],x[4]),mc,AmericanOption(x[5],K)), x)
@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1]),equitySpotData(x[2],x[3],x[4]),mc,BarrierOptionDownOut(x[5],K,D)), x)
@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1]),equitySpotData(x[2],x[3],x[4]),mc,AsianFloatingStrikeOption(x[5])), x)
@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1]),equitySpotData(x[2],x[3],x[4]),mc,AsianFixedStrikeOption(x[5],K)), x)