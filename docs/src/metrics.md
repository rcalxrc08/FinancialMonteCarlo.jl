# [Metrics](@id Metric)

The following type of *Metric* are supported from the package:

* `pricer`
* `delta`
* `rho`
* `variance`
* `confinter`

## Common Interface

Each Metric must implement its own *Metric* method; the general interface is the following:
```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::ZeroRate,mcConfig::MonteCarloConfiguration,abstractPayoff::FinancialMonteCarlo.AbstractPayoff)
variance(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::ZeroRate,mcConfig::MonteCarloConfiguration,abstractPayoff::FinancialMonteCarlo.AbstractPayoff)
confinter(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::ZeroRate,mcConfig::MonteCarloConfiguration,abstractPayoff::FinancialMonteCarlo.AbstractPayoff,alpha::Real=0.99)
delta(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::ZeroRate,mcConfig::MonteCarloConfiguration,abstractPayoff::FinancialMonteCarlo.AbstractPayoff,dS0::Real=1e-7)
rho(mcProcess::FinancialMonteCarlo.BaseProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoff::FinancialMonteCarlo.AbstractPayoff, dr::Real = 1e-7)
FinancialMonteCarlo.distribution(mcProcess::FinancialMonteCarlo.BaseProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, T::num_) where {num_ <: Number}
```