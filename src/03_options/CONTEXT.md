# Stage 03: payoffs, engines, smiles

The payoff hierarchy, the generic engines per family, the concrete payoffs and the smiles.

## Inputs

| Load | Why |
|:--|:--|
| `options.jl` | the aggregator |
| `options/payoff_base.jl` | the hierarchy, `@option_type`, `vanilla_payout`, `call_put_sign`, `maturity`, `ControlVariates`, sampling checks for payoffs |
| `options/engines.jl` | includes the engines |
| `options/engines/general_european_engine.jl` | `payoff` of European payoffs, `values_at_maturity` |
| `options/engines/general_american_engine.jl` | Longstaff–Schwartz: `exercise_columns`, `lsm_regressor`, `lsm_coefficients`, `solve_3x3` |
| `options/engines/general_path_dependent_engine.jl` | `payoff_by_rows`, `path_range`, `row_means` (barriers, Asians) |
| `options/payoffs.jl` | includes the concrete payoffs |
| `options/payoffs/spot.jl`, `options/payoffs/forward.jl`, `options/payoffs/european_option.jl`, `options/payoffs/binary_european_option.jl` | `Spot`, `Forward`, `EuropeanOption`, `BinaryEuropeanOption` |
| `options/payoffs/barrier_do_option.jl`, `options/payoffs/barrier_di_option.jl`, `options/payoffs/barrier_uo_option.jl`, `options/payoffs/barrier_ui_option.jl`, `options/payoffs/double_barrier_option.jl` | barriers: `is_alive` only |
| `options/payoffs/american_option.jl`, `options/payoffs/binary_american_option.jl`, `options/payoffs/bermudan_option.jl` | American, binary American, Bermudan |
| `options/payoffs/asian_fixed_strike_option.jl`, `options/payoffs/asian_floating_strike_option.jl` | Asians: `asian_payout` |
| `options/payoffs/basket/n_european_option.jl` | `EuropeanOptionND` |
| `options/smiles_base.jl` | `OptionSmile`, `smile_payoffs`, `smile_to_options`, `validate_abstract_smile` |
| `options/smiles.jl` | `EuropeanOptionSmile`, `BinaryEuropeanOptionSmile`, `payout_strike`, `adapt_structure` |
| `references/payoffs.md` | the hierarchy, engines, Longstaff–Schwartz, smiles, adding a payoff |
| `references/autodiff.md` | `path_range` inactive, no LAPACK, positions under Enzyme |

## Process

- A new payoff: struct + `payout`/`is_alive`/`asian_payout`; a payoff holding arrays needs `adapt_structure`; a new
  Asian one a place in the `Union` of the CUDA extension; a row in the payoff bestiary.
- A new smile family: `payout_strike` and its struct.
- Payoffs expose `.T`.
- `options/payoffs/basket/n_average_forward.jl` is an untracked experiment, not included.

## Outputs

`payoff(S, payoff, zeroCurve, mcConfig, T)` returns discounted per-path payoffs (a matrix for smiles). Checked by
`test/engines.jl`, `test/models/*.jl` (`check_prices`), `test/models/barrier_closed_form.jl`.

## Verify (Kaimon)

Setup: `using FinancialMonteCarlo, FinancialToolbox`; `bs = BlackScholesProcess(0.2, Underlying(100.0, 0.0))`,
`mc = MonteCarloConfiguration(100_000, 1; sampling = FinancialMonteCarlo.AntitheticMC())`.

| Check | Expect |
|:--|:--|
| `pricer(bs, ZeroRate(0.01), mc, EuropeanOption(1.0, 100.0)) - blsprice(100.0, 100.0, 0.01, 1.0, 0.2, 0.0)` | within a few standard errors (observed -0.050 with seed 0; absolute error below 0.15) |
| `length(pricer(bs, ZeroRate(0.01), mc, FinancialMonteCarlo.EuropeanOptionSmile(1.0, [90.0, 100.0], [true, true])))` | `2` |
