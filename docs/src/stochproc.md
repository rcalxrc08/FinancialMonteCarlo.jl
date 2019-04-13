# [Stochastic Processes](@id stochproc)

The package support the following processes:

* `Black Scholes Process`
* `Heston Process`
* `Kou Process`
* `Merton Process`
* `Variance Gamma Process`
* `Normal Inverse Gaussian Process`
* `Lognormal Mixture Process`

To be noticed that in order to work also with package DifferentialEquations.jl, the base class of the process is from it.

## Common Interface

A single method is implemented for each process, which provide the simulation output.

### Simulation

Each process behave in its own different way but returning the same kind of object after simulation,
the generic interface for simulating is the following:
```@docs
simulate(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::numb,monteCarloMode::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode()) where {numb<:Number}
```

