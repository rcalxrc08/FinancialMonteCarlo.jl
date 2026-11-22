# Stage 04: contract algebra and portfolios

_`Position`, `MarketDataSet`, `Portfolio`._ (Layer 3 reference.)

- `Position{P, W}` holds its legs as two vectors, `payoffs::Vector{P}` and `weights::Vector{W}`, in insertion order (a
  position has few legs, so a payoff is found by `isequal`; an `OrderedDict` there was untypable for Enzyme, #40):
  payoffs combine with `+`, `-`, `*`, `/` through `combine` (payoff type `typejoin`, weight `promote_type`, `foreach` over
  the tuple of terms); a weight that is exactly zero (`is_null_weight`: integers, floats, rationals only, never dual
  numbers) drops its leg. Payoffs compare by structure: `hash`/`isequal`/`==` on `Base.typename` and
  `parameter_values`, so `K = 100` and `K = 100.0` are one leg.
- A position is priced path by path: `payoff(S, position, ...)` is the weighted sum of its legs' discounted payoffs, so
  `pricer`, `variance`, `confinter`, `delta`, `rho` take positions; `maturity(position)` starts its fold from a zero of
  the maturity type (an `Int8` in the fold made a union that Enzyme rejects).
- `MarketDataSet` and `Portfolio` wrap `OrderedDict{UnderlyingKey, ...}`, keys a `String` or a tuple of names for
  multivariate models and baskets (`("ENI", "AAPL") → copula`); `"A_B"` is deprecated (`underlying_key`).
  `pricer(::MarketDataSet, ..., ::Portfolio)` seeds once and simulates the models in order from one stream
  (`@set mcConfig.reseed = false`), so the underlyings are independent (#38); `positions_on` gives each multivariate
  model its positions per subset of components (`compute_indices(n)`: the nonempty subsets of `1:n` by size, then
  lexicographically).
- Printing is `show(io, ::MIME"text/plain", ...)`, no `display` override; `capture_stdout` in the tests checks the
  overloads that write to stdout.
- A `Position`, a `Portfolio` or a `Vector{AbstractPayoff}` holds abstract types, so their prices are not inferred.
- Open: #65 (Enzyme on a Position with legs of different payoff types), #82 (arrays of payoffs and positions price a
  maturity off the time grid at the nearest node), #46 and #23 (Gaussian copula: cost of `quantile!`, path continuity).
- No multicurrency support.
