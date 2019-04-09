# FinancialMonteCarlo.jl Documentation

```@contents
```

## Functions

```@docs
payoff(S::AbstractMatrix{numtype},amPayoff::AmericanOption,spotData::equitySpotData,T1::num2=amPayoff.T) where{numtype,num2<:Number}
payoff(S::AbstractMatrix{num},asianFixedStrikePayoff::AsianFixedStrikeOption,spotData::equitySpotData,T1::num2=asianFixedStrikePayoff.T) where{num,num2<:Number}
payoff(S::AbstractMatrix{num},asianFloatingStrikePayoff::AsianFloatingStrikeOption,spotData::equitySpotData,T1::num2=asianFloatingStrikePayoff.T) where{num,num2<:Number}
payoff(S::AbstractMatrix{num},barrierPayoff::BarrierOptionDownIn,spotData::equitySpotData,T1::num2=barrierPayoff.T) where{num,num2<:Number}
payoff(S::AbstractMatrix{num},barrierPayoff::BarrierOptionDownOut,spotData::equitySpotData,T1::num2=barrierPayoff.T) where{num,num2<:Number}
payoff(S::AbstractMatrix{num},barrierPayoff::BarrierOptionUpIn,spotData::equitySpotData,T1::num2=barrierPayoff.T) where{num,num2<:Number}
payoff(S::AbstractMatrix{num},barrierPayoff::BarrierOptionUpOut,spotData::equitySpotData,T1::num2=barrierPayoff.T) where{num,num2<:Number}
payoff(S::AbstractMatrix{num},amPayoff::BinaryAmericanOption,spotData::equitySpotData,T1::num2=amPayoff.T) where{num,num2<:Number}
payoff(S::AbstractMatrix{num},euPayoff::BinaryEuropeanOption,spotData::equitySpotData,T1::num2=euPayoff.T) where{num,num2<:Number}
payoff(S::AbstractMatrix{num},doubleBarrierPayoff::DoubleBarrierOption,spotData::equitySpotData,T1::num2=doubleBarrierPayoff.T) where{num,num2<:Number}
payoff(S::AbstractMatrix{num},euPayoff::EuropeanOption,spotData::equitySpotData,T1::num2=euPayoff.T) where{num,num2<:Number}
payoff(S::AbstractMatrix{num},optionData::Forward,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}
payoff(S::AbstractMatrix{num},spotData::equitySpotData,phi::Function,T::num2) where{num,num2<:Number}
payoff(S::AbstractMatrix{num},optionData::Forward,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}
```

## Index

```@index
```