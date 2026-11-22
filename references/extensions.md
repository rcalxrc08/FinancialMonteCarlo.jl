# Extensions

_What each weak dependency adds and the rules for extensions._ (Layer 3 reference.)

- One extension per weak dependency: CUDA (GPU models, `ext/FinancialMonteCarloCUDAExt/`), DifferentialEquations (price
  directly from an SDE), DualNumbers, EnzymeCore, Mooncake, VectorizedRNG, Distributed. TaylorSeries needs no extension.
- Code there only exists when the user has loaded the corresponding package. An extension can only add methods, so
  types and functions it needs must be defined in `src/` (`CudaMode`, `MultiProcess`, `MonteCarloDiffEqModel`,
  `HostAdaptor`, `default_nbatches`, `run_on_batches`, `array_type`, `mode_adaptor`, `seed_device!`), and it imports the
  unexported names of the package explicitly.
- The CUDA extension changes behaviour only through the mode (`array_type`, `mode_adaptor`, `seed_device!`) and the array
  type (`CuArray`/`CuMatrix` methods), never through the configuration type. It holds `device_draw` (Gamma,
  InverseGaussian), the jump kernel, `draw_normals!`, the `Union` of Asian payoffs. CUDA is not a test dependency, so
  nothing runs that extension on CI (`test/cuda/`, locally).
- Mooncake: `set_seed!` zero-derivative primitive; `rrule!!`/`frule!!` of the kernels. EnzymeCore: inactive `set_seed!`
  and scalar draws; `augmented_primal`/`reverse` of the kernels (`autodiff.md`).
- Distributed: `MultiProcess` batches (`run_on_batches` with `pmap`).
- DifferentialEquations: `MonteCarloDiffEqModel`; #64 asks for DifferentialEquations 8.
