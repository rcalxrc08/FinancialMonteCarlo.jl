# Stage 03: payoffs and smiles

_The payoff hierarchy, the engines, Longstaff–Schwartz, smiles._ (Layer 3 reference; the catalogue is
`docs/src/manual/payoff_bestiary.md`.)

## Hierarchy (`payoff_base.jl`)

- `SingleNamePayoff` → `EuropeanPayoff`, `AmericanPayoff`, `BermudanPayoff`, `PathDependentPayoff` (→ `BarrierPayoff`,
  `AsianPayoff`); separately `BasketPayoff`.
- Generic engines in `options/engines/` implement `payoff(S, ...)` per family; a concrete payoff usually only defines
  its struct and a scalar `payout(ST, payoff)` (European) or path-based equivalent.
- Options made of positive numeric fields plus `isCall` are defined with `@option_type Name Family T K ...`, and
  `vanilla_payout` holds the call/put formula; barrier options only define `is_alive(S, payoff)`.
- Payoffs must expose `.T` (used by `maturity`).
- Barrier and Asian payoffs are row reductions (`payoff_by_rows`, from `barrier_payout(ST, lo, hi, payoff)` /
  `asian_payout(ST, A, payoff)`; `is_alive` takes the path or its `PathRange`), so a new barrier option needs nothing and
  a new Asian one its `asian_payout` and a place in the `Union` of the CUDA extension. `path_range` (the minimum and
  maximum for the barriers) is inactive under AD, and Asian means are a matrix-vector product (`row_means`).

## American and Bermudan

One Longstaff–Schwartz `payoff`, differing in `exercise_columns`. Its 3×3 normal equations are solved on scalars by
`solve_3x3` (Gaussian elimination with partial pivoting, as LAPACK): no allocations, and faster than `lu!` in every
backend, Enzyme included, which has no rule for LAPACK's `getrf`. The regression is on the return from the spot,
`(S - S0) / S0` (`lsm_regressor`), not on the price: on prices around 100 the columns 1, S, S² are almost collinear and
the fitted continuation values lost two digits (#35, tested against a `BigFloat` fit in `test/engines.jl`); the sums
and the elimination use `muladd`, the continuation value `evalpoly`.

## Smiles (`smiles_base.jl`, `smiles.jl`)

- `EuropeanOptionSmile`, `BinaryEuropeanOptionSmile`, aliases of `OptionSmile`, are not `AbstractPayoff`s (a smile has
  one value per strike, so it stays out of `Position`), but have the same interface: `payout(ST, smile)` returns one
  value per strike, `payoff(S, smile, ...)` a matrix with one row per path and one column per strike, and `pricer`,
  `variance`, `confinter`, `delta`, `rho` one value per strike, without building the options.
- A smile and its single option share one formula, `payout_strike(ST, K, isCall, OptionType)`, so a new smile family
  needs only that method and its struct.
- Strikes and flags are `AbstractVector`s, and `smile_payoffs` is one broadcast of the values at maturity against
  `transpose(K)` on any device (no `reshape`); one-dimensional inputs are typed `AbstractVector` throughout (Bermudan
  dates, curve nodes, `set_parameters!`).
- In `CudaMode` all the strikes run in one kernel and these metrics return a `CuArray` (`confinter` a `Vector` of
  tuples); `pricer` and `payoff_moments` move the smile to the mode first (`to_mode`), so a host smile is priced on the
  GPU and a GPU smile on the host without errors.
- `Adapt.adapt_structure` of a smile rebuilds it with the unchecked inner constructor (`Val(:unchecked)`): moving the
  arrays cannot invalidate it, and on the GPU `minimum(K)` was a device sync per metric call (#73).
- `smile_to_options` remains for what really needs the options (VibratoMonteCarlo's `singular_points`).
- #83 proposes a term structure of smiles for the pricer interface.

## Adding a payoff

A struct subtyping the right payoff family with its parameters (`T`, `K`, barriers…) and a `payoff`/`payout` method;
all metrics then work automatically. A payoff holding arrays needs `Adapt.adapt_structure`. A row in
`docs/src/manual/payoff_bestiary.md` with its Black–Scholes closed form (from FinancialToolbox).
