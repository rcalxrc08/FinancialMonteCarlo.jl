# FinancialMonteCarlo.jl Documentation

```@contents
```

## Functions

```@docs
payoff(S::AbstractMatrix{num},optionData::AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}
simulate(mcProcess::BaseProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::numb,monteCarloMode::MonteCarloMode=standard,parallelMode::BaseMode=SerialMode()) where {numb<:Number}
struct AmericanOption{num1,num2<:Number}<:AmericanPayoff
```

## Index

```@index
```