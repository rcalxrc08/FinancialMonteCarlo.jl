# Payoffs

A payoff struct holds the terms of a contract: maturity `T` and, depending on the contract,
strike `K`, barriers and whether it is a call (`isCall = true`, the default) or a put.

| Family | Payoffs |
|:--|:--|
| Linear | [`Spot`](@ref), [`Forward`](@ref) |
| European | [`EuropeanOption`](@ref), [`BinaryEuropeanOption`](@ref) |
| Barrier | [`BarrierOptionDownOut`](@ref), [`BarrierOptionDownIn`](@ref), [`BarrierOptionUpOut`](@ref), [`BarrierOptionUpIn`](@ref), [`DoubleBarrierOption`](@ref) |
| Asian | [`AsianFixedStrikeOption`](@ref), [`AsianFloatingStrikeOption`](@ref) |
| Early exercise | [`AmericanOption`](@ref), [`BinaryAmericanOption`](@ref), [`BermudanOption`](@ref) |
| Basket | [`EuropeanOptionND`](@ref) |
| Smiles | [`FinancialMonteCarlo.EuropeanOptionSmile`](@ref), [`FinancialMonteCarlo.BinaryEuropeanOptionSmile`](@ref) |

```julia
EuropeanOption(1.0, 100.0)                # 1y ATM call
EuropeanOption(1.0, 100.0, false)         # 1y ATM put
BarrierOptionDownOut(1.0, 100.0, 90.0)    # down and out call, barrier at 90
DoubleBarrierOption(1.0, 100.0, 80.0, 120.0)
BermudanOption([0.25, 0.5, 0.75, 1.0], 100.0)  # exercise dates
```

Constructors check that maturities, strikes and barriers are positive.

## How payoffs are evaluated

The generic [`payoff`](@ref) function maps the matrix of simulated paths to the vector of
discounted payouts, one per path. Each family has a generic *engine* implementing it, so a concrete
payoff only needs to describe its payout:

- **European** payoffs look at the column of the paths at maturity and call `payout(S_T, payoff)`;
- **path dependent** payoffs (barriers, Asians) receive the whole path of each simulation;
- **American** and **Bermudan** payoffs are priced with the Longstaff-Schwartz regression, using
  the exercise value `payout(S_t, payoff)` at every step (or exercise date).

Because the paths are shared, several payoffs with different maturities can be priced on the same
simulation: the paths go up to the largest maturity and each payoff picks the time steps it needs.
For this to be exact, maturities should fall on the time grid `0:T/Nstep:T`.

## Smiles

A smile groups European (or binary European) options with the same maturity and different strikes:

```julia
FinancialMonteCarlo.EuropeanOptionSmile(1.0, [90.0, 100.0, 110.0], [false, true, true])
```

`pricer` returns the vector of prices.

## Positions

Payoffs can be added, subtracted and multiplied by numbers, see
[Portfolios and multivariate](@ref).

See [Payoffs](@ref api-payoffs) in the API reference.
