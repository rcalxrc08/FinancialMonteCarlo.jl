using FinancialMonteCarlo, Test
@show "Test Parameters"
S0=100.0;
K=100.0;
r=0.02;
Tneg=-1.0;
d=0.01;

Nsim=10000;
Nstep=30;
sigma=0.2; 
McConfig=MonteCarloConfiguration(Nsim,Nstep);
toll=0.8;

spotData1=ZeroRateCurve(r);

@show "Test Geometric Brownian Motion Parameters"
drift=0.0
@test_throws(ErrorException,simulate(GeometricBrownianMotion(sigma,drift,Underlying(S0,d)),spotData1,McConfig,Tneg));
@test_throws(ErrorException,GeometricBrownianMotion(-sigma,drift,Underlying(S0,d)))