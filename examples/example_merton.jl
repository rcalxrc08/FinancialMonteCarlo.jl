using MonteCarlo
include("runner.jl")

Nsim=10000;
Nstep=30;
sigma=0.2;
lam=5.0; 
mu1=0.02; 
sigma1=0.05;
const dict=Dict{String,Number}("sigma"=>sigma, "lambda" => lam, "muJ" => mu1, "sigmaJ" => sigma1)
mc=MonteCarloBaseData(dict,Nsim,Nstep);
runnerMonteCarlo(MertonProcess(),mc)