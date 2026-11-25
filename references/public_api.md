# Public API

_The interface shared with the other packages and the main exports._ (Layer 3 reference; `docs/src/api.md` lists every
exported or documented name.)

## Interface shared with FinancialFFT, VibratoMonteCarlo, FinancialCalibrator

`pricer(model, underlying, zeroCurve, method, payoff)` with five arguments (since 0.10, D-014: the models hold only
their parameters; a multivariate model takes a tuple of underlyings; `pricer(market, zeroCurve, method, portfolio)`
for a `MarketDataSet`, which holds the underlyings by name), the mode inside the method (read with `execution_mode`,
changed with `with_mode`), device placement through `to_mode`/`mode_adaptor` and dispatch on `AbstractCudaMode`. A
change of these names is a breaking change for those packages too. VibratoMonteCarlo also reuses the jump helpers
(`JumpList`, `add_jump!`, `jump_quantile`, `on_sobol_paths`) and `smile_to_options`; FinancialCalibrator the curves,
`set_parameters!`/`get_parameters` and their Zygote rules.

## Exports (checked 2026-11-22 with `names(FinancialMonteCarlo)`)

| Group | Names |
|:--|:--|
| metrics | `pricer`, `variance`, `confinter`, `delta`, `rho`, `simulate`, `payoff`, `→` |
| configuration and market | `MonteCarloConfiguration`, `ZeroRate`, `Underlying` |
| engines | `BrownianMotion`, `GeometricBrownianMotion`, `GeometricBrownianMotionVec`, `SubordinatedBrownianMotion` |
| models | `BlackScholesProcess`, `BachelierProcess`, `HestonProcess`, `KouProcess`, `MertonProcess`, `VarianceGammaProcess`, `NormalInverseGaussianProcess`, `LogNormalMixture`, `ShiftedProcess`, `GaussianCopulaNVariateProcess`, `GaussianCopulaNVariateLogProcess`, `MonteCarloDiffEqModel` |
| payoffs | `Spot`, `Forward`, `EuropeanOption`, `BinaryEuropeanOption`, `BarrierOptionDownIn`, `BarrierOptionDownOut`, `BarrierOptionUpIn`, `BarrierOptionUpOut`, `DoubleBarrierOption`, `AmericanOption`, `BinaryAmericanOption`, `BermudanOption`, `AsianFixedStrikeOption`, `AsianFloatingStrikeOption`, `EuropeanOptionND` |
| portfolios | `MarketDataSet` |

Not exported, used qualified (`FinancialMonteCarlo.X`) by users and the other packages: the modes (`SerialMode`,
`MultiThreading`, `MultiProcess`, `CudaMode`, `AbstractCudaMode`), the samplings (`StandardMC`, `AntitheticMC`,
`SobolMC`, `ControlVariates`), `execution_mode`, `with_mode`, `to_mode`, `Curve`, `ImpliedCurve`, `underlying(market, key)`, the smiles
(`EuropeanOptionSmile`, `BinaryEuropeanOptionSmile`), `Position`, `Portfolio`, `distribution`, `get_parameters`,
`set_parameters!`, `integral`. `docs/src/api.md` lists every exported or documented name.
