using MonteCarlo
include("runner.jl")

Nsim=10000;
Nstep=30;
sigma=0.2;
lam=5.0; 
mu1=0.02; 
sigma1=0.05;
mc=MonteCarloBaseData(Nsim,Nstep);
runnerMonteCarlo(MertonProcess(sigma,lam,mu1,sigma1),mc)