using Test, FinancialMonteCarlo,Statistics;
@show "Black Scholes Model"
S0=100.0;
K=100.0;
r=[0.00,0.02];
T=1.0;
d=0.01;
r=0.02;
D=90.0;

Nsim=100000;
Nstep=30;
sigma=0.2;
mc=MonteCarloConfiguration(Nsim,Nstep);
toll=0.8

rfCurve=ZeroRate(r);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)
Model_base=BrownianMotion(sigma,r-d-sigma*sigma*0.5,Underlying(0.0,0.0))
S=S0.*exp.(simulate(Model_base,mc,T));

Payoff=exp(-FinancialMonteCarlo.integral(r,T))*S[:,end];
@show Price=mean(Payoff)
Model=BlackScholesProcess(sigma,Underlying(S0,d));
@show FwdPrice=pricer(Model,rfCurve,mc,FwdData);