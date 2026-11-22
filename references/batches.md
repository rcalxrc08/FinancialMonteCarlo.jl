# Stage 06: batched modes

_`MultiThreading`, `MultiProcess`._ (Layer 3 reference.)

- Batched modes split into serial batch configurations: `src/06_batches/batches.jl` holds all they share, and a mode only
  adds a `run_on_batches` method (`multi_threading.jl` for threads; `MultiProcess` in the Distributed extension).
- `batch_results` converts the results to the type inference gives (`fetch` and `pmap` return `Any`, which would make
  every batched metric type unstable).
- `pricer` averages the batches, `variance`/`confinter` pool their statistics, `simulate` joins their paths, so results
  do not depend on threads or workers.
- Batch `i` gets a copy of `rng` with seed `seed + i` (or `seed[i]`); each batch holds whole antithetic pairs
  (`batch_sizes`).
- Specializing the batch helpers on `f` (`batch_results`, `average_on_batches`, compiled for `f::Function`) was measured
  and gains nothing: the thread spawns dominate.
