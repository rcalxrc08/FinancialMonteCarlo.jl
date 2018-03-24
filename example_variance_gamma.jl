include("montecarlo.jl")
include("runner.jl")

Nsim=10000;
Nstep=30;
sigma=0.2;
theta1=-0.03; 
k1=0.16;
const dict=Dict{String,Number}("sigma"=>sigma, "theta" => theta1, "k" => k1)
mc=MonteCarloBaseData(dict,Nsim,Nstep);
runnerMonteCarlo(VarianceGammaProcess(),mc)