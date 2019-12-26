# [Metrics](@id Metric)

The following type of *Metric* are supported from the package:

* `European Options`
* `Asian Options`
* `Barrier Options`
* `American Options`

## Common Interface

Each Metric must implement its own *Metric* method; the general interface is the following:
```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::ZeroRateCurve,mcConfig::MonteCarloConfiguration,abstractPayoff::FinancialMonteCarlo.AbstractPayoff)
variance(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::ZeroRateCurve,mcConfig::MonteCarloConfiguration,abstractPayoff::FinancialMonteCarlo.AbstractPayoff)
confinter(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::ZeroRateCurve,mcConfig::MonteCarloConfiguration,abstractPayoff::FinancialMonteCarlo.AbstractPayoff,alpha::Real=0.99)
delta(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::ZeroRateCurve,mcConfig::MonteCarloConfiguration,abstractPayoff::FinancialMonteCarlo.AbstractPayoff,dS0::Real=1e-7)
```