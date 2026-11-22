# Dependencies

_Why each dependency exists, hard or weak, and its compat floor._ (Layer 3 reference.)

| Package | Hard/weak | Why | Compat |
|:--|:--|:--|:--|
| `Accessors` | hard | `@set` on immutable models and methods (`with_mode`, delta/rho bumps) | 0.1 |
| `Adapt` | hard | `to_mode`, `adapt_structure` of smiles | 3, 4 |
| `ChainRulesCore` | hard | `@ignore_derivatives`, the two rrules (#80) | 1 |
| `Distributions`, `StatsFuns` | hard | subordinator laws, quantiles | 0.25; 1, 2 |
| `MuladdMacro` | hard | `@muladd` in kernels | 0.2.7 |
| `OrderedCollections` | hard | `MarketDataSet`, `Portfolio` in insertion order | 1.6, 2 |
| `Sobol` | hard | `SobolMC` (`sobol_skip!` sets the state in closed form) | 1.0.2 |
| `StaticArrays` | hard | `SVector` curve nodes (#79, since 0.9) | 1.5 |
| `LinearAlgebra`, `Random`, `Statistics` | hard | stdlib | 1 |
| `CUDA` | weak | GPU models, device draws | 5, 6 |
| `StochasticDiffEq` | weak | `MonteCarloDiffEqModel` (DifferentialEquations 8 has no SDE solvers, #64) | 6, 7 |
| `Distributed` | weak | `MultiProcess` | 1 |
| `DualNumbers` | weak | dual-number methods | 0.6 |
| `EnzymeCore` | weak | inactive draws, kernel rules | 0.8 |
| `Mooncake` | weak | zero-derivative `set_seed!`, kernel rules | 0.5 |
| `VectorizedRNG` | weak | an alternative rng | 0.2 |

Julia floor 1.10 (CI 1.10 and 1.13.3).

Upstream workarounds (cite issue and version; recheck on every upstream release, remove when fixed):

| Workaround | Upstream issue | Affected versions | Last checked |
|:--|:--|:--|:--|
| CI pinned to Julia 1.13.3 | EnzymeAD/Enzyme.jl#4051 (fixed in Enzyme 0.13.220) | Enzyme < 0.13.220 on Julia 1.13.4 | 2026-11-18 |
| `set_runtime_activity` for positions on Julia 1.10 | Enzyme on 1.10 (#66) | Julia 1.10 | 2026-11-19 |
| `@test_broken` Enzyme/Mooncake in `CudaMode` | JuliaGPU/CUDA.jl#3332, EnzymeAD/Enzyme.jl#3741, #4094 (#68) | CUDA 6.5.1, Enzyme 0.13.221 | 2026-11-19 |
| `promote`/`convert` instead of `oftype`; benchmarks without splats | FluxML/Zygote.jl#1751, #1752, #1753, #1754 | Zygote 0.7.18 | 2026-11-19 |
