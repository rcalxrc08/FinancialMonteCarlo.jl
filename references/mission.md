# Mission

_Why the package exists and how to tell whether it succeeds._ (Layer 3 reference: stable, changes only with a decision
in `decisions.md`.)

## Purpose

Monte Carlo pricing and risk of equity derivatives, single-name and multivariate (Gaussian copula), with one interface
(`pricer(model, underlying, zeroCurve, method, payoff)`) that FinancialFFT, VibratoMonteCarlo and FinancialCalibrator extend with
their own methods; every model and payoff differentiable by every Julia AD back-end and runnable on CPU threads,
processes and CUDA.

## Goals

- A model catalogue (Black–Scholes, Bachelier, Heston, Kou, Merton, VG, NIG, lognormal mixture, shifted, copulas) and a
  payoff catalogue (European, binary, barriers, American, Bermudan, Asian, baskets, smiles, positions, portfolios).
- Samplings: standard, antithetic, Sobol, control variates; results independent of threads and workers.
- AD with rules on the hot kernels (Mooncake, Enzyme), Zygote on the structs (#80).

## Constraints

- Julia ≥ 1.10; registered in General.
- No Float64 hard-coding; one `simulate!` body per model; extensions only add methods.
- The five-argument `pricer` and the mode-in-method contract are shared with three packages.

## Success criteria

- The suite passes on CI (1.10, 1.13.3); `test/cuda/` locally; Aqua clean; JET finds only the Sobol barrier.
- Prices within the tolerances of the FinancialToolbox closed forms.

## Non-goals

- Simulations under Zygote (only the structs); multicurrency; Fourier pricing (FinancialFFT); calibration
  (FinancialCalibrator).
