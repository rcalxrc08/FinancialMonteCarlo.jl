include("montecarlo.jl")
include("runner.jl")

Nsim=10000;
Nstep=30;
sigma=0.2; 
const dict=Dict{String,Number}("sigma"=>sigma)
mc=MonteCarloBaseData(dict,Nsim,Nstep);
runnerMonteCarlo(BlackScholesProcess(),mc)