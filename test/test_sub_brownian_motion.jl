using FinancialMonteCarlo, Distributions, Test;
@show "Test Parameters"
S0=100.0;
r=0.02;
d=0.01;
T=1.0;
Nsim=10000;
Nstep=30;

sigma=0.2; 

McConfig=MonteCarloConfiguration(Nsim,Nstep);
spotData1=ZeroRateCurve(r);

@show "Test Subordinated Brownian Motion Parameters"
drift=0.0;

@test_throws(ErrorException,simulate(SubordinatedBrownianMotion(sigma,drift,InverseGaussian(1.0,1.0),Underlying(S0,d)),spotData1,McConfig,-T));
@test_throws(ErrorException,SubordinatedBrownianMotion(-sigma,drift,InverseGaussian(1.0,1.0),Underlying(S0,d)))