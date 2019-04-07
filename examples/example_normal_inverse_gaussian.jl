using FinancialMonteCarlo
include("runner.jl")

Nsim=10000;
Nstep=30;
sigma=0.2;
theta1=-0.03; 
k1=0.16; 
mc=MonteCarloConfiguration(Nsim,Nstep);
runnerMonteCarlo(NormalInverseGaussianProcess(sigma,theta1,k1),mc)