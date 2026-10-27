# Market data

## Underlying

The spot price and the dividend yield of an asset are stored in an [`Underlying`](@ref), which
is part of every single name model:

```julia
Underlying(100.0)          # no dividends
Underlying(100.0, 0.01)    # 1% continuous dividend yield
Underlying(100.0, FinancialMonteCarlo.Curve([0.01, 0.015, 0.02], 1.0))  # time dependent dividends
```

## Zero rates

Discounting and the risk neutral drift come from a zero rate curve:

```julia
ZeroRate(0.02)                                      # constant rate
ZeroRate([0.00, 0.01, 0.02], 1.0)                   # instantaneous rates on a uniform grid up to 1.0
FinancialMonteCarlo.ImpliedZeroRate([0.01, 0.015], [0.5, 1.0])  # from zero rates quoted at given maturities
```

## Curves

Time dependent quantities are represented as curves of *instantaneous* rates, linearly interpolated
between the nodes. Models never need the instantaneous value itself, only its integral between two
times, which is what is used both for the drift of each simulation step and for discounting.

- [`FinancialMonteCarlo.Curve`](@ref) builds the curve directly from instantaneous rates;
- [`FinancialMonteCarlo.ImpliedCurve`](@ref) builds the instantaneous curve that reproduces given
  *zero* rates, i.e. such that ``\int_0^{T_i} r(t)\,dt = r_i T_i``. This is what you want when you
  have market quoted zero rates or dividend yields.

There is no support for multi-currency.

See [Market data](@ref api-market-data) in the API reference for the full signatures.
