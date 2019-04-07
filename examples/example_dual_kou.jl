using FinancialMonteCarlo
include("runner.jl")

Nsim=10000;
Nstep=30;
sigma=0.2; 
p=0.3; 
lam=5.0; 
lamp=30.0; 
lamm=20.0; 
mc=MonteCarloConfiguration(Nsim,Nstep);
runnerMonteCarloDual(KouProcess(sigma,lam,p,lamp,lamm),mc)