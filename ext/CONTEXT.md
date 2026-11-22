# Extensions

One extension per weak dependency; the package loads and works without any of them.

## Inputs

| Load | Why |
|:--|:--|
| `FinancialMonteCarloCUDAExt/FinancialMonteCarloCUDAExt.jl` | `array_type`, `mode_adaptor`, `seed_device!` for `CudaMode`; `draw_normals!`, `device_draw` (Gamma, InverseGaussian), the jump kernel, Asian `Union` |
| `FinancialMonteCarloMooncakeExt.jl` | zero-derivative `set_seed!`; `rrule!!`/`frule!!` of the kernels |
| `FinancialMonteCarloEnzymeCoreExt.jl` | inactive `set_seed!` and scalar draws; `augmented_primal`/`reverse` of the kernels |
| `FinancialMonteCarloDualNumbersExt.jl` | dual-number methods |
| `FinancialMonteCarloStochasticDiffEqExt.jl` | simulation of `MonteCarloDiffEqModel` (SOSRI, `EnsembleThreads`) |
| `FinancialMonteCarloDistributedExt.jl` | `MultiProcess`: `run_on_batches` with `pmap` |
| `FinancialMonteCarloVectorizedRNGExt.jl` | VectorizedRNG generators |
| `references/extensions.md` | what lives in `src/` for the extensions, the CUDA rules |
| `references/autodiff.md` | the kernels and their rules, Enzyme and Mooncake constraints, GPU AD (#68) |
| `references/dependencies.md` | compat, upstream workarounds |

## Process

- An extension only adds methods; types and functions it needs are defined in `src/`, and it imports unexported names
  explicitly.
- CUDA changes behaviour only through the mode and the array type, never the configuration type; run `test/cuda/`
  locally (CI has no GPU).
- A new kernel's rules are added to both the Mooncake and the EnzymeCore extensions.

## Outputs

`Base.get_extension(FinancialMonteCarlo, :FinancialMonteCarlo<Dep>Ext) !== nothing` after `using <Dep>`. Checked by
`test/ad/mooncake.jl`, `test/ad/enzyme.jl`, `test/ad/dual_numbers.jl`, `test/models/diffeq.jl`,
`test/methods/multiprocess.jl`, `test/cuda/` (local).

## Verify (Kaimon)

| Check | Expect |
|:--|:--|
| `run_tests()` with the AD group | pass (Mooncake's first rule build ~2 min, Enzyme ~2.5 min) |
