using MonteCarlo
include("runner.jl")

Nsim=10000;
Nstep=30;
sigma=0.2; 
mc=MonteCarloBaseData(Nsim,Nstep);
runnerMonteCarloDual(BlackScholesProcess(sigma),mc)