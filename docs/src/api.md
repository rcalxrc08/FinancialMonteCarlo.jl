# API reference

```@contents
Pages = ["api.md"]
Depth = 2:2
```

## [Metrics](@id api-metrics)

```@docs
pricer
variance
confinter
delta
rho
FinancialMonteCarlo.distribution
```

## [Market data](@id api-market-data)

```@docs
Underlying
ZeroRate
ZeroRateCurve
FinancialMonteCarlo.ImpliedZeroRate
FinancialMonteCarlo.Curve
FinancialMonteCarlo.ImpliedCurve
```

## [Monte Carlo configuration](@id api-configuration)

```@docs
MonteCarloConfiguration
FinancialMonteCarlo.StandardMC
FinancialMonteCarlo.AntitheticMC
FinancialMonteCarlo.SobolMode
FinancialMonteCarlo.PrescribedMC
FinancialMonteCarlo.ControlVariates
FinancialMonteCarlo.SerialMode
FinancialMonteCarlo.MultiThreading
```

## [Models](@id api-models)

```@docs
simulate
BlackScholesProcess
HestonProcess
KouProcess
MertonProcess
VarianceGammaProcess
NormalInverseGaussianProcess
LogNormalMixture
ShiftedLogNormalMixture
GaussianCopulaNVariateProcess
GaussianCopulaNVariateLogProcess
```

### Engines

```@docs
BrownianMotion
FinancialMonteCarlo.BrownianMotionVec
GeometricBrownianMotion
GeometricBrownianMotionVec
SubordinatedBrownianMotion
FinancialMonteCarlo.SubordinatedBrownianMotionVec
```

### Parameters

```@docs
FinancialMonteCarlo.get_parameters
FinancialMonteCarlo.set_parameters!
```

## [Payoffs](@id api-payoffs)

```@docs
payoff
Spot
Forward
EuropeanOption
BinaryEuropeanOption
BarrierOptionDownOut
BarrierOptionDownIn
BarrierOptionUpOut
BarrierOptionUpIn
DoubleBarrierOption
AsianFixedStrikeOption
AsianFloatingStrikeOption
AmericanOption
BinaryAmericanOption
BermudanOption
EuropeanOptionND
FinancialMonteCarlo.EuropeanOptionSmile
FinancialMonteCarlo.BinaryEuropeanOptionSmile
```

## [Portfolios](@id api-portfolios)

```@docs
MarketDataSet
FinancialMonteCarlo.Position
FinancialMonteCarlo.Portfolio
→
```

## Index

```@index
Pages = ["api.md"]
```
