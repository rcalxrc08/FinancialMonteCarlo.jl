# FinancialMonteCarlo.jl Documentation

```@contents
```

## Functions

```@docs
payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}
simulate(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcBaseData::MonteCarloConfiguration,T::numb,monteCarloMode::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode()) where {numb<:Number}
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoff::FinancialMonteCarlo.AbstractPayoff,mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
struct AmericanOption{num1,num2<:Number}<:AmericanPayoff
```

## Index

```@index
```