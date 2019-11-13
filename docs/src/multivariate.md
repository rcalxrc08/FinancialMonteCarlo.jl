# [Payoffs](@id payoff)

The following type of *Payoff* are supported from the package:

* `European Options`
* `Asian Options`
* `Barrier Options`
* `American Options`

## Common Interface

Each payoff must implement its own *payoff* method; the general interface is the following:
```@docs
payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::ZeroRateCurve,T1::num2=optionData.T) where{num <: Number,num2 <: Number}
```
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
```