# [Multivariate Distributions](@id multivariates)

*Multivariate distributions* are the distributions whose variate forms are `Multivariate` (*i.e* each sample is a vector). Abstract types for multivariate distributions:

```julia
const MultivariateDistribution{S<:ValueSupport} = Distribution{Multivariate,S}

const DiscreteMultivariateDistribution   = Distribution{Multivariate, Discrete}
const ContinuousMultivariateDistribution = Distribution{Multivariate, Continuous}
```

## Common Interface

The methods listed as below are implemented for each multivariate distribution, which provides a consistent interface to work with multivariate distributions.

### Computation of statistics

```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
```

### Probability evaluation

```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
```
**Note:** For multivariate distributions, the pdf value is usually very small or large, and therefore direct evaluating the pdf may cause numerical problems. It is generally advisable to perform probability computation in log-scale.


### Sampling

```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
```

**Note:** In addition to these common methods, each multivariate distribution has its own special methods, as introduced below.


## Distributions

```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
```

## Addition Methods

### AbstractMvNormal

In addition to the methods listed in the common interface above, we also provide the following methods for all multivariate distributions under the base type `AbstractMvNormal`:

```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
```

### MvLogNormal

In addition to the methods listed in the common interface above, we also provide the following methods:

```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
```

It can be necessary to calculate the parameters of the lognormal (location vector and scale matrix) from a given covariance and mean, median or mode. To that end, the following functions are provided.

```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
```

## Internal Methods (for creating you own multivariate distribution)

```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
```
