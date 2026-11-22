# Stage 01: configuration and market data

Methods, modes and samplings; the configuration and its random state; rates, curves and underlyings; validation.

## Inputs

| Load | Why |
|:--|:--|
| `utils.jl` | the aggregator: includes the files below in order |
| `utils/checks.jl` | `check_domain`, `check_dimension`, `check_argument` |
| `utils/curves.jl` | `AbstractCurve`, `Curve{N, T1, T2}` on `SVector`s, `ImpliedCurve`, `CurveDifference`, `value_at`, `integral`, `incremental_integral` and the `rrule` of `integral` (`integral_partials`) |
| `utils/rng_number.jl` | `get_rng_type`: the number type of the draws |
| `utils/zero_rates.jl` | `ZeroRate`, `ZeroRateCurve`, `ImpliedZeroRate` |
| `utils/mc_mode.jl` | `AbstractMethod`, `AbstractMode`, `SerialMode`, `MultiProcess`, `CudaMode`, `AbstractCudaMode`, the samplings (`StandardMC`, `AntitheticMC`, `SobolMC`), `execution_mode`, `with_mode`, `mode_adaptor`, `to_mode`, `HostAdaptor`, `default_nbatches` |
| `utils/mc_config.jl` | `MonteCarloConfiguration`, `check_configuration` (`check_sampling`, `check_sampling_mode`, `check_reseed`, `check_mode`), `random_state`, `SeededRNG`, `host_rng`, `with_mode`, `SobolConfig`, `BatchedConfig` |
| `utils/underlying.jl` | `AbstractUnderlying`, `Underlying`, `UnderlyingScalar`, `UnderlyingVec` |
| `references/configuration.md` | methods, modes, samplings, data placement, curves |
| `references/conventions.md` | before adding a function |
| `references/autodiff.md` | isbits curves for Enzyme and the GPU; the Zygote rule of `integral` |

## Process

- A new compatibility restriction is a dispatch method of `check_configuration`, not a check in the metrics.
- Curves stay isbits (`SVector` nodes); between two curves only `+` and `-`.
- `with_mode` rebuilds `MonteCarloConfiguration`'s random state; never add a trailing `mode` argument anywhere.

## Outputs

Configurations, rates and underlyings used by every later stage. Checked by `test/structs.jl`, `test/methods/*.jl`,
`test/ad/chain_rules.jl` (`integral`), `test/ad/enzyme.jl` (curve underlyings).

## Verify (Kaimon)

Setup: `using FinancialMonteCarlo, Test`; `mc = MonteCarloConfiguration(1000, 10)`.

| Check | Expect |
|:--|:--|
| `FinancialMonteCarlo.execution_mode(mc)` | `SerialMode()` |
| `isbits(FinancialMonteCarlo.Curve([0.01, 0.02], [1.0, 2.0]))` | `true` |
| `MonteCarloConfiguration(1001, 10; sampling = FinancialMonteCarlo.AntitheticMC())` | an `ArgumentError` (odd `Nsim` with antithetic sampling) |
