```@raw html
<div class="fc-page-head">
<p class="fc-kicker"><span>Manual 07</span> Every payoff, one line each</p>
```

# Payoff bestiary

Every contract the package prices: what it pays, its constructor, and the engine that evaluates it. How payoffs are
evaluated, smiles and positions are explained in [Payoffs](@ref); this page is the catalogue.

```@raw html
<div class="fc-meta"><span>linear</span><span>European</span><span>barrier</span><span>Asian</span><span>early exercise</span><span>basket</span><span>smiles</span></div>
<nav class="fc-toc">
<a href="#Notation"><b>01</b>Notation</a>
<a href="#Linear"><b>02</b>Linear</a>
<a href="#European"><b>03</b>European</a>
<a href="#Barrier"><b>04</b>Barrier</a>
<a href="#Asian"><b>05</b>Asian</a>
<a href="#Early-exercise"><b>06</b>Early exercise</a>
<a href="#Basket"><b>07</b>Basket</a>
<a href="#Smiles"><b>08</b>Smiles</a>
<a href="#Combinations"><b>09</b>Combinations</a>
</nav>
</div>
```

## Notation

The paths are simulated on the grid ``t_k = k\,T / n``, ``k = 0, \dots, n``, with ``S_{t_0} = S_0``. A payoff of
maturity ``T`` pays at ``T`` and is discounted with ``D_T = e^{-\int_0^T r_s\,ds}``; ``ω = 1`` for a call
(`isCall = true`, the default) and ``ω = -1`` for a put, and ``(x)^+ = \max(x, 0)``. Barriers and averages are
monitored on the grid, so they converge to their continuous counterparts as the steps grow. Several payoffs priced
together share the paths, up to the longest maturity; maturities should fall on the grid.

| Payoff | Pays | Engine | Closed form under Black-Scholes |
|:--|:--|:--|:--|
| [`Spot`](@ref) | ``S_0`` today | none | ``S_0`` |
| [`Forward`](@ref) | ``S_T`` | European | ``S_0 e^{-\int_0^T d}`` |
| [`EuropeanOption`](@ref) | ``(ω (S_T - K))^+`` | European | `blsprice` |
| [`BinaryEuropeanOption`](@ref) | ``1_{ω (S_T - K) > 0}`` | European | `blsbin` |
| [`BarrierOptionDownOut`](@ref) | ``(ω (S_T - K))^+ 1_{\min_k S_{t_k} ≥ D}`` | path dependent | `blsdownout` (continuous monitoring) |
| [`BarrierOptionDownIn`](@ref) | ``(ω (S_T - K))^+ 1_{\min_k S_{t_k} < D}`` | path dependent | `blsdownin` (continuous monitoring) |
| [`BarrierOptionUpOut`](@ref) | ``(ω (S_T - K))^+ 1_{\max_k S_{t_k} ≤ U}`` | path dependent | `blsupout` (continuous monitoring) |
| [`BarrierOptionUpIn`](@ref) | ``(ω (S_T - K))^+ 1_{\max_k S_{t_k} > U}`` | path dependent | `blsupin` (continuous monitoring) |
| [`DoubleBarrierOption`](@ref) | ``(ω (S_T - K))^+ 1_{D ≤ \min_k S_{t_k},\ \max_k S_{t_k} ≤ U}`` | path dependent | |
| [`AsianFixedStrikeOption`](@ref) | ``(ω (A_T - K))^+`` | path dependent | geometric average only |
| [`AsianFloatingStrikeOption`](@ref) | ``(ω (S_T - A_T))^+`` | path dependent | |
| [`AmericanOption`](@ref) | ``(ω (S_τ - K))^+`` at the exercise ``τ ≤ T`` | Longstaff-Schwartz | |
| [`BinaryAmericanOption`](@ref) | ``1_{ω (S_τ - K) > 0}`` at the exercise ``τ ≤ T`` | Longstaff-Schwartz | |
| [`BermudanOption`](@ref) | ``(ω (S_τ - K))^+``, ``τ`` one of the dates | Longstaff-Schwartz | |
| [`EuropeanOptionND`](@ref) | ``(ω (\sum_i S^i_T - K))^+`` | basket | |

with ``A_T = \frac{1}{n + 1} \sum_{k=0}^n S_{t_k}`` the arithmetic average on the grid, the spot included.

## Linear

```julia
Spot()
Forward(T)
```

[`Spot`](@ref) is the underlying itself, worth ``S_0`` whatever the model and the method: in a position it is the
hedge leg (`EuropeanOption(T, K) - Δ * Spot()`). [`Forward`](@ref) pays ``S_T``, worth ``S_0 e^{-\int_0^T d}`` in every
model of the [Model bestiary](@ref): it checks the risk neutral drift of a model.

## European

```julia
EuropeanOption(T, K, isCall = true)
BinaryEuropeanOption(T, K, isCall = true)    # cash or nothing: pays 1
```

| Parameter | Meaning | Constraint |
|:--|:--|:--|
| `T` | maturity | ``T > 0`` |
| `K` | strike | ``K > 0`` |
| `isCall` | call (`true`) or put | |

Evaluated on the column of the paths at maturity only. Every model has a reference price for them (see the
**Closed form** column of the [Model bestiary](@ref)).

## Barrier

```julia
BarrierOptionDownOut(T, K, D, isCall = true)   # knocked out below D
BarrierOptionDownIn(T, K, D, isCall = true)    # knocked in below D
BarrierOptionUpOut(T, K, U, isCall = true)     # knocked out above U
BarrierOptionUpIn(T, K, U, isCall = true)      # knocked in above U
DoubleBarrierOption(T, K, D, U, isCall = true) # knocked out below D or above U
```

| Parameter | Meaning | Constraint |
|:--|:--|:--|
| `T`, `K`, `isCall` | as for the European option | |
| `D`, `U` (`barrier`) | lower and upper barrier | positive; ``D < U`` |

A European payout times whether the path is alive, from its minimum and maximum on the grid (spot and maturity
included). An in and the matching out option add up to the European one, path by path. Discrete monitoring
overprices a knock-out against the continuous barrier; the usual correction moves the barrier away from the spot
by a factor ``e^{0.5826\, σ \sqrt{T/n}}`` (Broadie, Glasserman and Kou): with it, the closed forms of FinancialToolbox
(`blsdownout`, `blsdownin`, `blsupout`, `blsupin`) match the Monte Carlo prices of the four options under Black-Scholes
within 2% (`test/models/barrier_closed_form.jl`). Their pathwise Greeks (any AD backend) are biased: the indicator
does not see the paths that cross the barrier when the spot or the volatility moves (#74); bump and reprice with a
common seed (`delta`) instead. On the GPU the minimum and maximum are row reductions,
so no path leaves the device.

## Asian

```julia
AsianFixedStrikeOption(T, K, isCall = true)    # pays (ω (A_T - K))⁺
AsianFloatingStrikeOption(T, isCall = true)    # pays (ω (S_T - A_T))⁺
```

The average ``A_T`` is arithmetic, on the ``n + 1`` points of the grid up to maturity, the spot included. On the
GPU it is a row reduction.

## Early exercise

```julia
AmericanOption(T, K, isCall = true)
BinaryAmericanOption(T, K, isCall = true)
BermudanOption(dates, K, isCall = true)        # dates: sorted vector, the last one the maturity
```

Exercisable at every step of the grid (American) or at the given dates (Bermudan, which should fall on the grid). The
exercise policy is estimated by the Longstaff-Schwartz regression of the continuation value on ``1, x, x^2``, with
``x = (S_t - S_0) / S_0``, over the paths in the money; the price is the mean discounted payout at the estimated exercise. A Bermudan with every step of the grid as a date is the
American on the same paths, and one with the maturity only is the European. American and Bermudan payoffs collect
the paths on the host, also in `CudaMode`.

## Basket

```julia
EuropeanOptionND(T, K, isCall = true)          # pays (ω (Σᵢ Sⁱ_T - K))⁺
```

On a multivariate model ([`GaussianCopulaNVariateProcess`](@ref), [`GaussianCopulaNVariateLogProcess`](@ref)),
whose components it sums at maturity; see [Portfolios and multivariate](@ref).

## Smiles

```julia
FinancialMonteCarlo.EuropeanOptionSmile(T, K, isCall)         # K, isCall: vectors
FinancialMonteCarlo.BinaryEuropeanOptionSmile(T, K, isCall)
```

The European options of one maturity and several strikes, priced on the same paths: the metrics return one value
per strike, and in `CudaMode` all the strikes run in one kernel. See [Payoffs](@ref).

## Combinations

Payoffs on the same underlying add, subtract and scale into a [`FinancialMonteCarlo.Position`](@ref) (a straddle is
`EuropeanOption(T, K) + EuropeanOption(T, K, false)`), and positions on several underlyings form a
[`FinancialMonteCarlo.Portfolio`](@ref): see [Portfolios and multivariate](@ref). A new payoff needs only its
struct and its payout, see [Extending the package](@ref).
