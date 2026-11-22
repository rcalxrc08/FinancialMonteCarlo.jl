# Workflows: commands, Kaimon, CI, releases, docs deployment

_How to run, test, benchmark, document and release._ (Layer 3 reference.)

## Kaimon session

- `start_session(project_path = "C:/Users/Nicola/.julia/dev/FinancialMonteCarlo")`; `kaimon.toml` pins Julia 1.13.
- First call: `ex(q = false, e = "(VERSION, Base.active_project(), pathof(FinancialMonteCarlo))")`.
- AD back-ends are test dependencies: check them through `run_tests()` or a group in a separate process.
- Do not run `julia --project` scripts that precompile this package while the session is open (Revise keeps the old
  code); restart the session afterwards.

## Commands

```bash
julia --project=. -e 'using Pkg; Pkg.test()'                                   # full suite (coverage = true for .cov files)
julia --project=. -e 'using Pkg; Pkg.test(test_args = ["models", "ad"])'       # groups of TEST_GROUPS, case-insensitive
julia --project=test test/runtests.jl models
julia --project=test -e 'using Pkg; Pkg.instantiate()'                         # once
julia --project=test test/models/black_scholes.jl                              # one file
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'                    # once
julia --project=benchmark -t auto benchmark/run.jl                             # all groups
BENCHMARK_SECONDS=0.3 julia --project=benchmark -t auto benchmark/run.jl payoffs greeks
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl                                              # Documenter + Literate
julia -e 'using JuliaFormatter; format(".")'                                   # .JuliaFormatter.toml, margin=1200
```

## Remotes and CI

- The main remote is GitHub (`github.com/rcalxrc08/FinancialMonteCarlo.jl`, the `GitHub` remote); GitLab (`origin`,
  whose push URLs also include GitHub) is a legacy mirror, not pushed since 0.5.2: push to `GitHub`.
- GitHub Actions: CI (Julia 1.10 and 1.13.3), Benchmarks, Documentation, CompatHelper, TagBot. Since 2026-11-14 CI pins
  1.13.3 instead of 1: Enzyme aborted on Julia 1.13.4 (EnzymeAD/Enzyme.jl#4051, fixed in Enzyme 0.13.220, which
  FinancialToolbox already relies on): put '1' back once the Enzyme tests pass on 1.13.4 (`state/STATUS.md`).
- The CI and Documentation runs of 54ecd0b (2026-11-11) stayed `queued` for hours with Actions reported operational,
  and `gh run cancel` answered that they were completed: the next push started fresh runs normally.

## Release

1. Bump `version` in `Project.toml` on master (breaking changes: minor bump while 0.x), wait for CI to be green.
2. Register by commenting on the release commit on GitHub, with release notes (sections: Breaking changes, New features,
   Performance, Bug fixes): `gh api repos/rcalxrc08/FinancialMonteCarlo.jl/commits/<sha>/comments -F body=@notes.md`,
   where `notes.md` starts with `@JuliaRegistrator register` then `Release notes:`. AutoMerge merges in ~20 minutes.
3. **Tag by hand**, TagBot usually does not: `git tag -a vX.Y.Z -m "FinancialMonteCarlo vX.Y.Z" <sha>`,
   `git push GitHub vX.Y.Z`, then `gh release create vX.Y.Z -R rcalxrc08/FinancialMonteCarlo.jl --title vX.Y.Z
   --notes-file <notes without the registrator lines> --verify-tag`.
4. Check that `https://rcalxrc08.github.io/FinancialMonteCarlo.jl/stable/` serves the new version.

## Docs deployment

- GitHub Pages builds from GitHub Actions (`build_type: workflow`); the `source.branch: gh-pages` in
  `gh api repos/rcalxrc08/FinancialMonteCarlo.jl/pages` is a leftover. There is no `gh-pages` branch (#52): the `docs`
  job of `Documentation.yml` builds the manual of master into `docs/build` (the site root), checks out the highest `v*`
  tag (`git ls-remote --tags | sort -V`) into `release/` and builds its manual into `docs/build/stable` (with
  `GITHUB_REPOSITORY` unset, so the `deploydocs` of the `make.jl` of releases up to 0.7 does nothing), adds
  `docs/redirect/404.html`, and uploads `docs/build` as the Pages artifact, deployed by `actions/deploy-pages`.
  `docs/make.jl` has no `deploydocs`, and its `canonical` is the site root.
- Old URLs: `dev/...` and `vX.Y.Z/...` no longer exist; `404.html` redirects them to the root (`dev/`) or `stable/` (a
  version folder). `stable/...` keeps working. Docs of older releases are built locally from their tag.
- A tag push runs only the `publish-tag` job, which starts the workflow on master (Pages refuses deployments from a tag;
  a `workflow_dispatch` started with `GITHUB_TOKEN` does start a run). Pull requests build the manual of master only.

## Sub-environments

`test/`, `docs/` and `benchmark/` load the package through `[sources]`; each lists the package itself in `[compat]`
(bump on breaking releases). FinancialToolbox is a dependency of `docs/` and `test/`.
