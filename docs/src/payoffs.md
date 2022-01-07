# [Payoffs](@id payoff)

The following type of *Payoff* are supported from the package:

* `European Options`
* `Asian Options`
* `Barrier Options`
* `American Options`
* `Bermudan Options`

## Common Interface

Each payoff must implement its own *payout* method; that will be called by the following general interface:
```@docs
payoff(S::AbstractMatrix, p::FinancialMonteCarlo.AbstractPayoff, r::FinancialMonteCarlo.AbstractZeroRateCurve, mc::FinancialMonteCarlo.AbstractMonteCarloConfiguration, T1::Number = p.T)
```
Depending on the type of AbstractPayoff, a different type of S must be specified. Look into the engines for more details.
The following option structs are provided:
```@docs
BinaryEuropeanOption
AsianFloatingStrikeOption
BarrierOptionDownOut
AsianFixedStrikeOption
DoubleBarrierOption
BarrierOptionUpOut
BarrierOptionUpIn
BinaryAmericanOption
EuropeanOption
Forward
BarrierOptionDownIn
AmericanOption
Spot
EuropeanOptionND
FinancialMonteCarlo.BermudanOption
```