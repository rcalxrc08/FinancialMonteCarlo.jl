# Truncated Distributions

The package provides the `Truncated` type to represented truncated distributions, which is defined as below:

```julia
struct Truncated{D<:UnivariateDistribution,S<:ValueSupport} <: Distribution{Univariate,S}
    untruncated::D      # the original distribution (untruncated)
    lower::Float64      # lower bound
    upper::Float64      # upper bound
    lcdf::Float64       # cdf of lower bound
    ucdf::Float64       # cdf of upper bound

    tp::Float64         # the probability of the truncated part, i.e. ucdf - lcdf
    logtp::Float64      # log(tp), i.e. log(ucdf - lcdf)
end
```

A truncated distribution can be constructed using the constructor `Truncated` as follows:

```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
```

Many functions, including those for the evaluation of pdf and sampling, are defined for all truncated univariate distributions:

- [`payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}`](@ref)
- [`payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}`](@ref)
- [`payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}`](@ref)
- [`payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}`](@ref)
- [`payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}`](@ref)
- [`payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}`](@ref)
- [`payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}`](@ref)
- [`payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}`](@ref)
- [`payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}`](@ref)
- [`payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}`](@ref)
- [`payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}`](@ref)
- [`payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}`](@ref)
- [`payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}`](@ref)
- [`payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}`](@ref)
- [`payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}`](@ref)
- [`payoff(S::AbstractMatrix{num},optionData::FinancialMonteCarlo.AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}`](@ref)

Functions to compute statistics, such as `mean`, `mode`, `var`, `std`, and `entropy`, are not available for generic truncated distributions. Generally, there are no easy ways to compute such quantities due to the complications incurred by truncation.
However, these methods are supported for truncated normal distributions.

```@docs
pricer(mcProcess::FinancialMonteCarlo.BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{FinancialMonteCarlo.AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::FinancialMonteCarlo.BaseMode=SerialMode())
```
