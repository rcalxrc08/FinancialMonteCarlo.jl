# FinancialMonteCarlo.jl Documentation

```@contents
```

## Functions

```@docs
payoff(::AbstractMatrix{numtype},::AmericanOption,::equitySpotData,::num2) where{numtype,num2<:Number}
payoff(S::AbstractMatrix{num},euPayoff::EuropeanOption,spotData::equitySpotData,T1::num2=euPayoff.T) where{num,num2<:Number}
```

## Index

```@index
```