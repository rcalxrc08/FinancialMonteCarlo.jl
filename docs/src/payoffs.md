# [Payoffs](@id payoff)

The following type of *Payoff* are supported from the package:

* `European Options`
* `Asian Options`
* `Barrier Options`
* `American Options`
* `Bermudan Options`

## Common Interface

Each payoff must implement its own *payout* method; the general interface is the following:
```@docs
payout(S,optionData::FinancialMonteCarlo.AbstractPayoff,spotData::FinancialMonteCarlo.AbstractZeroRateCurve,T1::num2=optionData.T) where{num <: Number,num2 <: Number}
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
```