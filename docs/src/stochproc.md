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
the generic interface for simulating is the following:
```@docs
simulate!(X,mcProcess::FinancialMonteCarlo.BaseProcess,spotData::ZeroRate,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: FinancialMonteCarlo.AbstractMonteCarloMethod, type4 <: FinancialMonteCarlo.BaseMode}
```
or
```@docs
simulate!(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::ZeroRate,mcBaseData::MonteCarloConfiguration{type1,type2,type3,type4},T::numb) where {numb <: Number, type1 <: Number, type2<: Number, type3 <: FinancialMonteCarlo.AbstractMonteCarloMethod, type4 <: FinancialMonteCarlo.BaseMode}
```
The following process are already implemented:
```@docs
BlackScholesProcess
BrownianMotion
GeometricBrownianMotion
HestonProcess
KouProcess
LogNormalMixture
MertonProcess
NormalInverseGaussianProcess
SubordinatedBrownianMotion
ShiftedLogNormalMixture
VarianceGammaProcess
```