# [Matrix-variate Distributions](@id matrix-variates)

*Matrix-variate distributions* are the distributions whose variate forms are `Matrixvariate` (*i.e* each sample is a matrix). Abstract types for matrix-variate distributions:

## Common Interface

Both distributions implement the same set of methods:

```@docs
payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}
```

## Distributions

```@docs
Wishart
InverseWishart
```

## Internal Methods (for creating your own matrix-variate distributions)

```@docs
payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}
```
