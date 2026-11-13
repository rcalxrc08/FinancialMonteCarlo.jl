```@raw html
<div class="fc-page-head">
<p class="fc-kicker"><span>Reference</span> Every documented name</p>
```

# API reference


```@raw html
<nav class="fc-toc">
<a href="#api-metrics"><b>01</b>Metrics</a>
<a href="#api-market-data"><b>02</b>Market data</a>
<a href="#api-configuration"><b>03</b>Monte Carlo configuration</a>
<a href="#api-models"><b>04</b>Models</a>
<a href="#api-payoffs"><b>05</b>Payoffs</a>
<a href="#api-portfolios"><b>06</b>Portfolios</a>
<a href="#Index"><b>07</b>Index</a>
</nav>
</div>
```
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
FinancialMonteCarlo.AbstractMethod
FinancialMonteCarlo.AbstractMonteCarloConfiguration
FinancialMonteCarlo.execution_mode
FinancialMonteCarlo.with_mode
FinancialMonteCarlo.AbstractSampling
FinancialMonteCarlo.StandardMC
FinancialMonteCarlo.AntitheticMC
FinancialMonteCarlo.SobolMC
FinancialMonteCarlo.ControlVariates
FinancialMonteCarlo.AbstractMode
FinancialMonteCarlo.SerialMode
FinancialMonteCarlo.MultiThreading
FinancialMonteCarlo.MultiProcess
FinancialMonteCarlo.CudaMode
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

### Jumps from uniforms

What `AntitheticMC` and `SobolMC` use to draw the jumps of `KouProcess` and `MertonProcess` from a
fixed number of uniforms per path (also used by VibratoMonteCarlo.jl).

```@docs
FinancialMonteCarlo.JumpUniforms
FinancialMonteCarlo.jump_count
FinancialMonteCarlo.jump_quantile
```

### Parameters

```@docs
FinancialMonteCarlo.get_parameters
FinancialMonteCarlo.set_parameters!
FinancialMonteCarlo.price_type
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
