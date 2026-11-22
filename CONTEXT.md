# CONTEXT.md: task routing (Layer 1)

Find the task, open the contract it points to, and load only the inputs listed there.

| Task | Contract |
|:--|:--|
| Configuration, modes, samplings, seeds, rates and curves, underlyings, validation helpers | `src/01_utils/CONTEXT.md` |
| A model, an engine, the simulation hooks, Sobol, jumps, the AD kernels and their adjoints | `src/02_models/CONTEXT.md` |
| A payoff, an engine (European, Longstaff–Schwartz, path-dependent), smiles | `src/03_options/CONTEXT.md` |
| `get_parameters`/`set_parameters!`, positions, portfolios, market data sets | `src/04_operations/CONTEXT.md` |
| `pricer`, `variance`, `confinter`, Greeks, control variates, output types | `src/05_metrics/CONTEXT.md` |
| Multi-threading and batched modes | `src/06_batches/CONTEXT.md` |
| CUDA, Mooncake, Enzyme, DualNumbers, StochasticDiffEq, Distributed, VectorizedRNG | `ext/CONTEXT.md` |
| Adding or running tests, groups, reference prices, the CUDA tests | `test/CONTEXT.md` |
| Documentation pages, bestiaries, tutorials, figures | `docs/CONTEXT.md` |
| The benchmark suite | `benchmark/CONTEXT.md` |
| The logo | `etc/CONTEXT.md` |
| Commands, CI, remotes, releases, docs deployment | `references/workflows.md` |
| What is being worked on, known issues | `state/STATUS.md` |
| Whether a change fits the package's goals; why a rule exists | `references/mission.md`, `references/decisions.md` |

Shared references (Layer 3), each loaded only when a contract lists it: `public_api.md`, `conventions.md`,
`autodiff.md`, `extensions.md`, `dependencies.md`, `workflows.md`, `testing.md`, `performance.md`, `docs_design.md`,
and per stage `configuration.md`, `models.md`, `payoffs.md`, `portfolios.md`, `metrics.md`, `batches.md`.

The stages communicate through Julia values: 01 builds the method and the market (`MonteCarloConfiguration`, rates,
underlyings); 02 returns the path matrix `Nsim × (Nstep+1)` from `simulate`; 03 returns discounted per-path payoffs
from `payoff(S, …)`; 04 combines payoffs into positions and portfolios; 05 reduces them to prices and statistics; 06
splits a batched method into serial batches and joins their results.
