# Getting Started

## Installation

The FinancialMonteCarlo package is available through the Julia package system by running `Pkg.add("https://gitlab.com/rcalxrc08/FinancialMonteCarlo.jl.git")`.
Throughout, we assume that you have installed the package.

## Hello World: European Option in Black Scholes Model

The following is the simplest example of use of FinancialMonteCarlo.jl.
```julia
#Import the Package
using FinancialMonteCarlo;
#Define Spot Datas
S0=100.0;
K=100.0;
r=0.02;
T=1.0;
d=0.01;
D=90.0;
#Define FinancialMonteCarlo Parameters
Nsim=10000;
Nstep=30;
#Define Model Parameters
σ=0.2;
#Build the Structs
mcConfig=MonteCarloConfiguration(Nsim,Nstep);
rfCurve=ZeroRate(r);

#Define The Options
EU_payoff=EuropeanOption(T,K)
#Define the Model
Model=BlackScholesProcess(σ,Underlying(S0,d));
#Price
EuPrice=pricer(Model,rfCurve,mcConfig,EU_payoff)
```

## Using Other Models and Options

The package contains a large number of models of three main types:

* `Ito Process`
* `Jump-Diffusion Process`
* `Infinity Activity Process`

And the following types of options:

* `European Options`
* `Asian Options`
* `Barrier Options`
* `American Options`

The usage is analogous to the above example.