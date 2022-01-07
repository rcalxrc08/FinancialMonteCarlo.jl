# [Stochastic Processes](@id stochproc)

The package support the following processes:

* `Black Scholes Process`
* `Heston Process`
* `Kou Process`
* `Merton Process`
* `Variance Gamma Process`
* `Normal Inverse Gaussian Process`
* `Lognormal Mixture Process`
* `Shifted Lognormal Mixture Process`

## Common Interface

A single method is implemented for each process, which provide the simulation output.

### Simulation

Each process behave in its own different way but returning the same kind of object after simulation,
the generic interfaces for simulating are the following:
```@docs
simulate(mcProcess::FinancialMonteCarlo.AbstractMonteCarloProcess,spotData::FinancialMonteCarlo.AbstractZeroRateCurve,mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, T::Number)
```
```@docs
simulate(mcProcess::FinancialMonteCarlo.AbstractMonteCarloEngine,mcBaseData::FinancialMonteCarlo.AbstractMonteCarloConfiguration, T::Number)
```
The following process are already implemented:
```@docs
BlackScholesProcess
BrownianMotion
FinancialMonteCarlo.BrownianMotionVec
GeometricBrownianMotion
FinancialMonteCarlo.GeometricBrownianMotionVec
HestonProcess
KouProcess
LogNormalMixture
MertonProcess
NormalInverseGaussianProcess
SubordinatedBrownianMotion
FinancialMonteCarlo.SubordinatedBrownianMotionVec
ShiftedLogNormalMixture
VarianceGammaProcess
```