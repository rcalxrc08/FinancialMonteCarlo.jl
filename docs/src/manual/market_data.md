```@raw html
<div class="fc-page-head">
<p class="fc-kicker"><span>Manual 02</span> Underlyings and curves</p>
```

# Market data


```@raw html
<div class="fc-meta"><span>Underlying</span><span>ZeroRate</span><span>Curve</span></div>
<nav class="fc-toc">
<a href="#Underlying"><b>01</b>Underlying</a>
<a href="#Zero-rates"><b>02</b>Zero rates</a>
<a href="#Curves"><b>03</b>Curves</a>
</nav>
</div>
```
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

Both take two vectors or two tuples. The nodes are stored as static vectors (StaticArrays), so a curve,
an underlying holding one and a model built on it are `isbits` values: Enzyme differentiates them
without runtime activity and a GPU kernel can capture them. The number of nodes is part of the type,
so a curve built from vectors has a type known only at run time: build it once, outside the code that
is differentiated or run in a loop, or from tuples, which give an inferred constructor.

Curves combine with numbers through `+`, `-`, `*` and `/`, applied to the values. Two curves can only
be added or subtracted: the result is a [`FinancialMonteCarlo.CurveDifference`](@ref), whose integral
is the difference of the integrals of the two curves, exact and without new nodes. This is how the
drift ``r - d`` is formed when both the rate and the dividend are curves. Curves are not multiplied or
divided by each other.

There is no support for multi-currency.

See [Market data](@ref api-market-data) in the API reference for the full signatures.
