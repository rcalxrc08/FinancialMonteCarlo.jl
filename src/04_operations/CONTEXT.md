# Stage 04: parameters and contract algebra

`get_parameters`/`set_parameters!`, display, positions, market data sets and portfolios.

## Inputs

| Load | Why |
|:--|:--|
| `io_utils.jl` | `get_parameters`, `set_parameters!` and its `rrule` (`parameter_tangent`), `parameter_fields`, `flat_parameter`, `parameter_vector`, `display`/`print` |
| `operations.jl` | the aggregator of the files below |
| `operations/models.jl` | `MarketDataSet`, `UnderlyingKey`, `underlying_key`, `underlying_names`, `number_type` |
| `operations/options.jl` | `Position`, `PositionTerm`, `combine`, `is_null_weight`, `maturity`, `payoff` of a position |
| `operations/portfolios.jl` | `Portfolio`, `positions_on` |
| `references/portfolios.md` | positions, portfolios, market data sets |
| `references/models.md` | `get_parameters`/`set_parameters!` rules ("Parameters") |
| `references/autodiff.md` | the Zygote rule of `set_parameters!`, Enzyme and positions |

## Process

- Parameters through the generated `parameter_names`/`parameter_values`, their number from `parameter_count(T)`; never
  a runtime Symbol.
- A model with non-scalar parameters (LogNormalMixture, ShiftedProcess) needs its own methods and a
  `parameter_tangent_values` method.
- Weights that are exactly zero drop their leg (never for dual numbers).

## Outputs

Parameter vectors, rebuilt models, positions and portfolios priced by stage 05. Checked by `test/models/parameters.jl`,
`test/ad/chain_rules.jl`, `test/portfolio/*.jl`.

## Verify (Kaimon)

Setup: `using FinancialMonteCarlo`; `bs = BlackScholesProcess(0.2, Underlying(100.0, 0.0))`.

| Check | Expect |
|:--|:--|
| `FinancialMonteCarlo.get_parameters(bs)` | `[0.2]` |
| `@inferred FinancialMonteCarlo.get_parameters(bs)` | a `Vector{Float64}` |
