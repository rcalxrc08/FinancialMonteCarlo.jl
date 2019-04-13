# [Univariate Distributions](@id univariates)

*Univariate distributions* are the distributions whose variate forms are `Univariate` (*i.e* each sample is a scalar). Abstract types for univariate distributions:

```julia
const UnivariateDistribution{S<:ValueSupport} = Distribution{Univariate,S}

const DiscreteUnivariateDistribution   = Distribution{Univariate, Discrete}
const ContinuousUnivariateDistribution = Distribution{Univariate, Continuous}
```

## Common Interface

A series of methods are implemented for each univariate distribution, which provide
useful functionalities such as moment computation, pdf evaluation, and sampling
(*i.e.* random number generation).

### Parameter Retrieval

**Note:** `params` are defined for all univariate distributions, while other parameter
retrieval methods are only defined for those distributions for which these parameters make sense.
See below for details.

```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
```

For distributions for which success and failure have a meaning,
the following methods are defined:
```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
```


### Computation of statistics

```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
```

### Probability Evaluation

### Sampling (Random number generation)
```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
```

## Continuous Distributions

```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
```

### Vectorized evaluation

Vectorized computation and inplace vectorized computation have been deprecated.
