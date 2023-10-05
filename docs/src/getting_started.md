```@meta
EditURL = "../../examples/getting_started.jl"
```

# Getting Started

## Installation
The FinancialMonteCarlo package is available through the Julia package system
by running `]add FinancialMonteCarlo`.
Throughout, we assume that you have installed the package.

## Basic syntax
The basic syntax for FinancialMonteCarlo is simple. Here is an example of pricing a European Option with a Black Scholes model:

````@example getting_started
using FinancialMonteCarlo
S0 = 100.0;
K = 100.0;
r = 0.02;
T = 1.0;
d = 0.01;
D = 90.0;
nothing #hide
````

Define FinancialMonteCarlo Parameters:

````@example getting_started
Nsim = 10000;
Nstep = 30;
nothing #hide
````

Define Model Parameters.

````@example getting_started
σ = 0.2;
nothing #hide
````

Build the MonteCarlo configuration and the zero rate.

````@example getting_started
mcConfig = MonteCarloConfiguration(Nsim, Nstep);
rfCurve = ZeroRate(r);
nothing #hide
````

Define The Option

````@example getting_started
EuOption = EuropeanOption(T, K)
````

Define the Model of the Underlying

````@example getting_started
Model = BlackScholesProcess(σ, Underlying(S0, d));
nothing #hide
````

Call the pricer metric

````@example getting_started
@show EuPrice = pricer(Model, rfCurve, mcConfig, EuOption);
nothing #hide
````

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

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

