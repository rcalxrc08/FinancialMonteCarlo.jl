# Conventions

_Coding standards that hold in every stage._ (Layer 3 reference.)

- Structs fully parametric with a promoted `numtype`; output types from `price_type`; no Float64 literals or buffers
  (`autodiff.md`).
- One `simulate!` body per model on the hooks of `paths.jl`; Sobol is the only separate body (`models.md`).
- Inner loops as kernels with closed-form adjoints and Mooncake/Enzyme rules; `@muladd` in kernels, Sobol loops, curves
  and device samplers.
- Validation by `check_domain`/`check_dimension`/`check_argument`, never `@assert` (#71); compatibility checks of the
  configuration in `check_configuration`, as dispatch methods.
- Every method carries its execution mode; never a trailing `mode` argument; inputs moved once with `to_mode`.
- Scalar draws through `draw_uniform`/`draw_exponential`/`draw_normal`, never `rand(host_rng(...))`; `foreach` over
  `path_sides`; no LAPACK in differentiated code.
- `get_parameters`/`set_parameters!` through the generated `parameter_names`/`parameter_values`, never runtime Symbols.
- Results inferred, also for arrays of payoffs of different types.
- Reuse before writing: closed forms from FinancialToolbox, never by hand.
- Each type file does its own `export`.
- Formatting: JuliaFormatter, `.JuliaFormatter.toml` (margin=1200, long one-line signatures expected).
