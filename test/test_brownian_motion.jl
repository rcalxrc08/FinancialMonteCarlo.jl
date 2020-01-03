using FinancialMonteCarlo, Test;
@show "Test Parameters"
S0=100.0;
K=100.0;
r=0.02;
Tneg=-1.0;
d=0.01;
T=-Tneg;
Nsim=10000;
Nstep=30;
mc1=MonteCarloConfiguration(Nsim,Nstep,FinancialMonteCarlo.AntitheticMC());
rfCurve2=FinancialMonteCarlo.ZeroRateCurve2([0.00,0.02],T);
sigma=0.2; 

McConfig=MonteCarloConfiguration(Nsim,Nstep);
rfCurve=ZeroRateCurve(r);

@show "Test Brownian Motion Parameters"
drift=0.0
@test_throws(ErrorException,BrownianMotion(-sigma,drift,Underlying(S0,d)))
@test_throws(ErrorException,BrownianMotion(-sigma,rfCurve2.r,Underlying(S0,d)))
@test_throws(ErrorException,simulate(BrownianMotion(sigma,drift,Underlying(S0,d)),McConfig,Tneg));
@test_throws(ErrorException,simulate(BrownianMotion(sigma,rfCurve2.r,Underlying(S0,d)),McConfig,Tneg));
@test_throws(ErrorException,simulate(BrownianMotion(sigma,drift,Underlying(S0,d)),mc1,Tneg));
@test_throws(ErrorException,simulate(BrownianMotion(sigma,rfCurve2.r,Underlying(S0,d)),mc1,Tneg));