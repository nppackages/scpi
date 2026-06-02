# Changelog

This file tracks notable changes made during the 2026 modernization and
refactoring effort for the `scpi` packages in R, Python, and Stata. Earlier
package history predates this changelog.

## 4.0.0 - In Progress

### Repository Modernization

- Reorganized the top-level README as a shared landing page for the R, Python,
  and Stata implementations.
- Updated package metadata to version `4.0.0`.
- Updated authorship metadata across platforms, with Matias D. Cattaneo as the
  maintainer where required by package infrastructure.
- Standardized references to the 2021 *Journal of the American Statistical
  Association* article, the 2025 *Journal of Statistical Software* article, and
  the 2027 *Review of Economics and Statistics* article.
- Added GitHub-oriented repository setup, including package release workflow
  groundwork and PyPI submission support.
- Aligned repository license text with the companion package setup.

### Local Workflow and Distribution

- Added local-only agent and workflow notes in `AGENTS.md`.
- Recorded local R, Python, and Stata 16 tool paths and package workflow notes.
- Established Stata 16 as the compatibility target for Stata compilation and
  tests.
- Updated the Stata distribution workflow to use a single `stata/lscpi.mlib`
  Mata library rather than distributing individual `.mo` files.
- Kept local auxiliary scripts under `script/`, which remains ignored by git.
- Added and refined local benchmark/check scripts for R, Python, and Stata
  numerical comparisons.
- Added an explicit local compatibility benchmark lane for running local source
  against older public Python dependency environments.
- Documented the Stata 16 embedded-Python setup: use base Python 3.10 and
  provide local package/dependency paths through `PYTHONPATH`.

### Stata

- Added `precision(single|double)` to Stata commands where generated-variable
  storage precision matters.
- Made `precision(double)` the default while preserving prior behavior through
  `precision(single)`.
- Updated Stata help files for the new precision option.
- Switched Stata wrappers to use direct SFI data access where possible.
- Removed a CSV round trip from Stata multi-unit plotting by passing plot data
  through Python/Stata data objects directly.
- Verified the Stata 16 smoke path with `scdatamulti`, `scpi`, and
  `scplotmulti`.

### R

- Precomputed repeated constants in uncertainty loops.
- Added safeguards around small simulation runs in diagonal parallel
  uncertainty code.
- Parallelized diagonal in-sample uncertainty simulation where appropriate.
- Tuned the diagonal parallel threshold to avoid overhead on small jobs.
- Batched multi-unit data stacking in `scdataMulti()` to avoid repeated object
  growth.
- Batched multi-unit pre-treatment fit assembly in `scest()`.
- Batched multi-unit inference design assembly in `scpi()`.
- Fixed time-effect aggregate interval row labels.
- Fixed sparse-matrix inference handling by using sparse-safe transposes.
- Fixed a unit-effect feature-design conformability issue when precomputed
  differenced post-treatment predictors are used with constants.
- Added regression tests for the unit-effect feature-design issue.
- Replaced formula-based weighted least squares inside the internal
  `shrinkage.EST()` ridge hyperparameter routine with an `lm.wfit()` helper
  that preserves the same weighted residual variance calculation with less
  setup overhead.
- Added an R reuse path allowing `scpi()` to accept a precomputed `scest`
  object and skip repeated point estimation. `scest` return objects now retain
  `Y.donors` in their data block, matching the documented return shape and
  making the object self-contained for inference reuse.

### Python

- Precomputed repeated constants in uncertainty loops.
- Batched multi-unit data concatenation in `scdataMulti()` to avoid repeated
  `pandas.concat()` growth.
- Fixed the Python counterpart of the unit-effect feature-design conformability
  issue.
- Fixed Python inference setup for Stata-fed multi-unit objects where
  precomputed differenced post-treatment predictors are present for only a
  subset of treated units.
- Fixed Python simultaneous prediction intervals when both columns of
  `w_bounds` are supplied, avoiding an uninitialized simulated-bounds object
  when in-sample uncertainty simulation is intentionally skipped.
- Added a regression test for the unit-effect feature-design issue.
- Added a Python reuse path allowing `scpi()` to accept precomputed `scest`
  output for single- and multi-unit objects and skip repeated point estimation.

### Numerical Checks and Testing

- Added local public-vs-local benchmark scaffolding for numerical drift checks.
- Compared selected old-vs-new R and Python scenarios for exact or near-exact
  numerical agreement after refactors.
- Verified R package tests with `testthat`.
- Verified R package checks with `R CMD check R\scpi --no-manual
  --no-build-vignettes`; current checks pass with the expected source-prep NOTE
  and restricted-network repository warnings.
- Verified Python changes with direct executable smoke tests in the local
  Python 3.10 benchmark environment. `pytest` is not currently installed in that
  environment.
- Verified the Stata 16 smoke workflow using base Python 3.10 plus local
  `PYTHONPATH` dependency paths.
- Ran local `single-smoke`, `single-illustration`, `multi-illustration`, and
  `restat-waveall` benchmark scenarios across R, Python, and Stata.
- Verified RESTAT WaveAll deterministic R/Python outputs within `5e-5`
  tolerance and RESTAT WaveAll Stata/Python outputs within `1e-6` tolerance.
- Observed public-vs-local inference drift in multi-unit interval/Sigma objects;
  deterministic RESTAT cross-platform fitted objects remain within tolerance.
- Verified that `origin/main` source and current local Python source produce
  identical RESTAT WaveAll snapshots under the same public Python dependency
  stack (`1e-8` tolerance over common numeric leaves).
- Verified that `origin/main` source and current local Python source produce
  identical multi-illustration snapshots under the same public Python dependency
  stack (`1e-8` tolerance over common numeric leaves).
- Verified that `origin/main` source and current local R source produce
  identical multi-illustration snapshots (`1e-8` tolerance over common numeric
  leaves).
- Confirmed that the RESTAT WaveAll R source-baseline run hits the aggregate
  time-effect interval dimname issue fixed locally.
- Isolated the larger RESTAT WaveAll Python inference drift to the installed
  PyPI package snapshot versus the GitHub source baseline, not to the
  modernization commits.
- Bounded RESTAT WaveAll Python/Stata dependency-stack drift at about `2e-6` in
  interval outputs when comparing public-compatible dependencies with the
  current local dependency stack.
- Verified that `origin/main` source and current local Stata source produce
  RESTAT WaveAll snapshots within `2e-6` tolerance under the same public
  Python dependency stack.
- Verified that Stata multi-illustration deterministic objects are stable
  against the `origin/main` source baseline, while interval matrices drift in
  that scenario; RESTAT remains the cleaner cross-platform numerical check.
- Verified the Stata 16 `scplotmulti` precision option with local Python
  dependencies: `precision(double)` stores generated plot variables as double,
  `precision(single)` preserves earlier single-precision generated storage, and
  invalid precision values return the expected error.
- Recorded preliminary single-rep speed checks: Python RESTAT WaveAll improved
  from about `14.50s` to `9.93s` against the GitHub source baseline under the
  same public dependency stack, and Python multi-illustration improved from
  about `2.35s` to `2.17s`.
- Recorded preliminary single-rep Stata speed checks: RESTAT WaveAll improved
  from about `22.99s` to `16.03s` against the GitHub source baseline under the
  same public Python dependency stack, and to about `14.40s` with the current
  local dependency stack.
- Ran 3-rep timing benchmarks for RESTAT WaveAll and multi-illustration. Median
  Python RESTAT runtime improved from `14.263s` (`origin/main`) to `10.186s`
  with current source under public-compatible dependencies (`1.40x`) and
  `8.771s` with current dependencies (`1.63x`).
- In the same 3-rep suite, median Python multi-illustration runtime improved
  from `2.454s` to `2.351s` under public-compatible dependencies (`1.04x`) and
  `2.244s` with current dependencies (`1.09x`).
- Median Stata RESTAT runtime improved from `19.632s` (`origin/main`) to
  `14.676s` with current source under public-compatible dependencies (`1.34x`)
  and `13.114s` with current dependencies (`1.50x`).
- Median Stata multi-illustration runtime improved modestly from `4.916s` to
  `4.704s` under public-compatible dependencies (`1.05x`) and `4.720s` with
  current dependencies (`1.04x`), with the same interval-matrix drift noted
  above.
- Median R multi-illustration runtime was `1.600s` for `origin/main` and
  `1.740s` locally in the 3-rep suite, while numerical outputs remained
  identical at `1e-8`; this small scenario is not yet evidence of an R speedup.
- The R RESTAT `origin/main` source baseline still fails with the known
  aggregate time-effect dimname issue, while current local R RESTAT completed
  with median runtime `5.890s`.
- Profiled current local R on multi-illustration and RESTAT WaveAll scenarios.
  The remaining R hotspots are repeated `b.est()` calls inside `scest()` and
  diagonal in-sample uncertainty simulation: in the RESTAT profile with
  `sims=10`, `scest()` took about `3.90s`, repeated `b.est()` calls accounted
  for about `3.70s`, and `insampleUncertaintyGetDiag()` took about `1.88s`.
- Tested a direct `ECOSolveR` point-estimation fast path for simplex problems.
  Loose solver tolerances caused unacceptable numerical drift; tight tolerances
  restored old-path equality at `1e-8` on single-smoke, multi-illustration, and
  RESTAT WaveAll checks, but did not improve runtime, so the experiment was not
  retained.
- Tested using `OSQP` for the lasso subproblems inside `shrinkage.EST()` and
  rejected it because RESTAT `Q2` values drifted materially and runtime did not
  improve.
- Verified the `lm.wfit()` shrinkage helper against the previous local R path:
  RESTAT WaveAll (`sims=10`) and multi-illustration (`sims=20`) matched at
  `1e-8` over common numeric leaves. Single-run total time moved from about
  `8.07s` to `8.03s` on RESTAT and from about `3.76s` to `3.65s` on
  multi-illustration.
- Verified the R `scest` reuse path with deterministic bounds in single- and
  multi-unit tests; reused and direct paths matched at `1e-8`. In a small
  deterministic timing check, median bounded `scpi()` time moved from about
  `0.03s` to `0.00s` in the single-unit case and from about `0.07s` to
  `0.02s` in the multi-unit case after reusing an existing `scest` object.
- Verified R tests with `devtools::test('R/scpi')`: 138 passing tests, no
  failures or warnings. Verified `R CMD check R\scpi --no-manual
  --no-build-vignettes` with the expected source-prep NOTE and restricted
  network repository warnings.
- Verified the Python `scest` reuse path with deterministic-bounds single- and
  multi-unit tests; direct and reused paths matched at `1e-8`. In a small
  5-rep Germany timing check with `w_bounds` and `sims=10`, mean single-unit
  `scpi()` time moved from about `0.044s` to `0.037s` after reusing an
  existing `scest` object (`1.19x`).
- Profiled current local Python on multi-illustration (`sims=20`) and RESTAT
  WaveAll (`sims=10`) scenarios. In multi-illustration, `scpi()` took about
  `2.36s`, with ECOS simulation solves accounting for about `1.50s` of self
  time across 2,052 small solver calls. In RESTAT WaveAll, `scdataMulti()` took
  about `3.27s` and `scpi()` about `8.97s`; remaining hotspots include ECOS
  simulation solves, out-of-sample quantile regression, and pseudo-inverse/SVD
  work.
- Tested storing Python quantile-regression predictions as numeric `float`
  arrays instead of object arrays. Multi-illustration and RESTAT WaveAll
  snapshots matched at `1e-8`, but single-run timing did not improve, so the
  experiment was not retained.

### Known Follow-Ups

- Run the broader public-vs-local numerical drift benchmark suite across R,
  Python, and Stata.
- Continue profiling repeated R point-estimation calls; any replacement for
  CVXR setup in `b.est()` needs strict old-path numerical validation and a real
  speed gain before it is kept.
- Continue profiling remaining inference hotspots after drift checks are in
  place.
- Install or provision `pytest` in the Python benchmark environment if full
  Python test-suite execution is desired locally.
- Continue compiling Stata help files to PDF before distribution.
