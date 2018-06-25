using MonteCarlo
include("runner.jl")

Nsim=10000;
Nstep=30;
sigma=0.2; 
mc=MonteCarloConfiguration(Nsim,Nstep);
runnerMonteCarlo(BlackScholesProcess(sigma),mc)