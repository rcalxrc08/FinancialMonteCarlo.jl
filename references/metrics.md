# Stage 05: metrics

_The pricing pipeline, control variates, variance and intervals, output types._ (Layer 3 reference.)

## Core pipeline (`metrics/pricer.jl`)

Every metric reduces to: `pricer(process, zeroCurve, mcConfig, payoff)` → `set_seed!` → `simulate(process, zeroCurve,
mcConfig, maturity(payoff))` returns an `Nsim × (Nstep+1)` path matrix → `payoff(S, payoff, zeroCurve, mcConfig, T)`
returns discounted per-path payoffs → `mean`. `pricer` is overloaded for: a single payoff, an array of payoffs, a
`Position`, an `AbstractSmile`, a vectorial process with `Array{Position}`, and a `MarketDataSet` with a `Portfolio`.
`delta`, `rho`, `variance`, `confinter`, `distribution` build on top (delta/rho are bump-and-reprice using
`Accessors.@set` on the immutable model structs).

## Control variates (`pricer_cv.jl`)

`ControlVariates(variate, pilot, sampling; price)` takes the exact discounted price of the control (`control_price`: the
given one, the forward's closed form, else the pilot mean, which cancels most of the reduction, #75); the pilot only
estimates the coefficient `c`.

## Variance and intervals (`variance.jl`, `confinterval.jl`)

The moments are those of the independent samples (`independent_samples`, `sample_count`): the paths, or with
`AntitheticMC` the averages of the pairs, path i and i + Nsim/2 (#76: the pairs are correlated, positively for an even
payoff, where treating them as independent gave intervals 30% too narrow). `variance` is per path (`per_path_variance`:
the sample variance of the samples times `paths_per_draw`), so `variance / Nsim` is the variance of the price in every
sampling; `confinter` uses `sample_count - 1` degrees of freedom; the batches pool the moments of their samples (each
batch holds whole pairs, `batch_sizes`). With `SobolMC` neither is an error estimate (docstrings).

## Output types (`output_type.jl`)

Output element types come from `price_type(inputs...)`, a number type (`typeof(zero(T1) + zero(T2) + …)` over the
inputs' `number_type`s), used to allocate path/price arrays of the right (possibly Dual) type: never hardcode `Float64`
buffers. It works on types only (no array traversal, so GPU arrays and ReverseDiff tapes are untouched);
`BaseProcess{T}` and `AbstractPayoff{T}` carry their type, anything else needs a `number_type` method (there is no
catch-all, so a missing one is a `MethodError`). `get_matrix_type(mcConfig, process, T)` takes that type.

Results must be inferred also for an array of payoffs of different types (abstract element type): take the element type
from `price_type` and annotate the comprehension (`pricer`, `payoff_moments`).

## Known biases

- #82: arrays of payoffs and positions price a maturity off the time grid at the nearest node (5-7% error).
- #74: pathwise (AD) Greeks of barrier options are biased (`@test_broken` in `test/models/barrier_closed_form.jl`).
- #55: Heston Monte Carlo is biased by ~10% when the Feller condition is badly violated.
- First- and second-order AD "works" for every payoff, including ones where it is meaningless (binary, barrier);
  second-order results are mostly useless.
- #28: European payoffs do not need the whole path matrix (performance).
