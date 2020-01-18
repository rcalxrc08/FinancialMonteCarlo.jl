using FinancialMonteCarlo, Test
@show "Test Parameters"
S0=100.0;
K=100.0;
r=0.02;
Tneg=-1.0;
T=-Tneg;
d=0.01;

Nsim=10000;
Nstep=30;
sigma=0.2; 
McConfig=MonteCarloConfiguration(Nsim,Nstep);
toll=0.8;
rfCurve2=ZeroRate([0.00,0.02],T);
rfCurve=ZeroRate(r);

@show "Test Geometric Brownian Motion Parameters"
drift=0.0
@test_throws(ErrorException,GeometricBrownianMotion(-sigma,drift,Underlying(S0,d)))
@test_throws(ErrorException,GeometricBrownianMotion(-sigma,rfCurve2.r,Underlying(S0,d)))

