using FinancialMonteCarlo, Test;
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
spotData1=equitySpotData(r,d);

@show "Test Brownian Motion Parameters"
drift=0.0
@test_throws(ErrorException,BrownianMotion(-sigma,drift,S0))
@test_throws(ErrorException,simulate(BrownianMotion(sigma,drift,S0),spotData1,McConfig,Tneg));