using MonteCarlo, Base.Test;
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
spotData1=equitySpotData(S0,r,d);

@show "Test Subordinated Brownian Motion Parameters"
drift=0.0;
dt=T/Nstep;
sub=dt*ones(Nsim,Nstep)

@test_throws(ErrorException,simulate(SubordinatedBrownianMotion(sigma,drift),spotData1,McConfig,-T,sub));
@test_throws(ErrorException,SubordinatedBrownianMotion(-sigma,drift))


@test_throws(ErrorException,SubordinatedBrownianMotion(-sigma,drift))
subw=dt*ones(1,1);
@test_throws(ErrorException,simulate(SubordinatedBrownianMotion(sigma,drift),spotData1,McConfig,T,subw));