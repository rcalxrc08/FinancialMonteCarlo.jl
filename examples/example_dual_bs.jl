using FinancialMonteCarlo
include("runner.jl")

Nsim=10000;
Nstep=30;
sigma=0.2; 
mc=MonteCarloConfiguration(Nsim,Nstep);
runnerMonteCarloDual(BlackScholesProcess(sigma),mc)