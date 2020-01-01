using Test, FinancialMonteCarlo,Statistics;
@show "Black Scholes Model"
S0=100.0;
K=100.0;
r=[0.00,0.02];
T=1.0;
d=0.01;
D=90.0;

Nsim=100000;
Nstep=30;
sigma=0.2;
mc=MonteCarloConfiguration(Nsim,Nstep);
toll=0.8

rfCurve=FinancialMonteCarlo.ZeroRateCurve2(r,T);

FwdData=Forward(T)
EUData=EuropeanOption(T,K)

S=S0.*exp.(simulate(BrownianMotion(sigma,rfCurve.r.-(sigma*sigma*0.5),Underlying(0.0,0.0)),rfCurve,mc,T));

Payoff=exp(-FinancialMonteCarlo.integral(rfCurve,T))*S[:,end];
Price=mean(Payoff)