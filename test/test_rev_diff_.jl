using FinancialMonteCarlo,ReverseDiff;
@show "Black Scholes Model"
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;

Nsim=100;
Nstep=30;
sigma=0.2; 
p=0.3; 
lam=5.0; 
lamp=30.0; 
lamm=20.0;
mc=MonteCarloConfiguration(Nsim,Nstep);
mc1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.AntitheticMC());
toll=0.8;

rfCurve=ZeroRate(r);

Model=BlackScholesProcess(sigma,Underlying(S0,d));

x=Float64[sigma,S0,r,d,T]

@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1],Underlying(x[2],x[4])),ZeroRate(x[3]),mc,Forward(x[5])), x)
@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1],Underlying(x[2],x[4])),ZeroRate(x[3]),mc,EuropeanOption(x[5],K)), x)
@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1],Underlying(x[2],x[4])),ZeroRate(x[3]),mc,AmericanOption(x[5],K)), x)
@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1],Underlying(x[2],x[4])),ZeroRate(x[3]),mc,BarrierOptionDownOut(x[5],K,D)), x)
@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1],Underlying(x[2],x[4])),ZeroRate(x[3]),mc,AsianFloatingStrikeOption(x[5])), x)
@show ReverseDiff.gradient(x->pricer(BlackScholesProcess(x[1],Underlying(x[2],x[4])),ZeroRate(x[3]),mc,AsianFixedStrikeOption(x[5],K)), x)