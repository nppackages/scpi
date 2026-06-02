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
- Removed a duplicate `Qtools::rrq()` fit from the internal R quantile
  regression prediction helper used by `e.method = "qreg"` and `e.method =
  "all"`.
- Extended `scdataMulti()` to accept unit-specific `features`, `cov.adj`,
  `constant`, and `cointegrated.data` specifications, matching the flexibility
  already available in the Stata and Python interfaces.
- Fixed R multi-unit `P.diff` handling for mixed cointegration settings so
  units that do not require precomputed differenced predictors are represented
  as `NULL` during inference instead of causing a column-name mismatch.
- Added regression coverage for unit-specific `scdataMulti()` options.

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
- Pre-generated serial Python in-sample uncertainty simulation draws in
  `scpi_in()` while preserving the existing random stream.
- Aligned the Python multi-unit illustration script with the Stata example:
  shared options for the unit-time effect and unit-specific options for the
  unit and time aggregated effects.

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
- Tested avoiding dense `Omega` multiplication in R `u.sigma.est()`. A
  row-weighted `crossprod()` rewrite caused interval drift above `1e-8`; an
  order-preserving sparse-diagonal rewrite matched at `1e-8` on
  multi-illustration and RESTAT WaveAll but did not improve timing
  consistently, so neither experiment was retained.
- Verified the R `scest` reuse path with deterministic bounds in single- and
  multi-unit tests; reused and direct paths matched at `1e-8`. In a small
  deterministic timing check, median bounded `scpi()` time moved from about
  `0.03s` to `0.00s` in the single-unit case and from about `0.07s` to
  `0.02s` in the multi-unit case after reusing an existing `scest` object.
- Added R regression coverage for unbounded simulated inference through the
  `scpi(scest_object, ...)` reuse path. With the same RNG stream as the direct
  `scpi(scdata_object, ...)` workflow, multi-illustration matched at `1e-8`;
  local timing checks moved the `scpi()` segment from `2.23s` to `1.97s` on
  multi-illustration and from `6.93s` to `3.19s` on RESTAT WaveAll, with total
  runtime changes of about `1.04x` and `1.01x`, respectively, after accounting
  for separately timed `scest()`.
- Verified R tests with `devtools::test('R/scpi')`: 141 passing tests, no
  failures or warnings. Verified `R CMD check R\scpi --no-manual
  --no-build-vignettes` with the expected source-prep NOTE and restricted
  network repository warnings.
- Tested preallocating cross-unit simulation matrices in R
  `insampleUncertaintyGetDiag()` instead of growing them with repeated
  `cbind()` calls. Multi-illustration and RESTAT WaveAll snapshots matched at
  `1e-8`, but single-run timing showed no improvement, so the experiment was
  not retained.
- Replaced the double-fit R `Qtools::rrq()` helper path with a single fit while
  preserving warning handling. A focused old-vs-new helper check matched
  exactly (`max_abs = 0`); helper timing over 100 repetitions moved from about
  `0.95s` to `0.43s` (`2.21x`), and a package-level `e.method = "qreg"` smoke
  path completed successfully.
- Re-verified R tests with `devtools::test('R/scpi')`: 141 passing tests, no
  failures or warnings.
- Verified the Python `scest` reuse path with deterministic-bounds single- and
  multi-unit tests; direct and reused paths matched at `1e-8`. In a small
  5-rep Germany timing check with `w_bounds` and `sims=10`, mean single-unit
  `scpi()` time moved from about `0.044s` to `0.037s` after reusing an
  existing `scest` object (`1.19x`).
- Added Python regression coverage for unbounded simulated inference through
  the `scpi(scest_output, ...)` reuse path. With the same RNG stream as the
  direct `scpi(scdata_output, ...)` workflow, the multi-unit check matched at
  `1e-8`; local timing checks moved total runtime from about `2.20s` to
  `2.13s` on multi-illustration (`1.03x`) and from about `8.34s` to `8.23s`
  on RESTAT WaveAll (`1.01x`), with the RESTAT `scpi()` segment itself moving
  from about `6.79s` to `6.17s` (`1.10x`) after separately timing `scest()`.
- Profiled current local Python on multi-illustration (`sims=20`) and RESTAT
  WaveAll (`sims=10`) scenarios. In multi-illustration, `scpi()` took about
  `2.36s`, with ECOS simulation solves accounting for about `1.50s` of self
  time across 2,052 small solver calls. In RESTAT WaveAll, `scdataMulti()` took
  about `3.27s` and `scpi()` about `8.97s`; remaining hotspots include ECOS
  simulation solves, out-of-sample quantile regression, and pseudo-inverse/SVD
  work.
- Rewrote Python `u_sigma_est()` to compute `Sigma` through `VZ` and diagonal
  residual weights instead of multiplying through a dense `Omega` matrix, while
  preserving the existing dense `u_var` output. Multi-illustration and RESTAT
  WaveAll snapshots matched at `1e-8`; the `scpi()` segment improved modestly
  from `1.935s` to `1.921s` on multi-illustration and from `6.823s` to
  `6.790s` on RESTAT in single-run checks.
- Pre-generated normal draws for serial Python `scpi_in()` simulations.
  Multi-illustration and RESTAT WaveAll snapshots matched at `1e-8`; in
  single-run checks, the `scpi()` segment moved from about `2.79s` to `2.76s`
  on multi-illustration and from about `9.54s` to `9.42s` on RESTAT WaveAll.
- Checked the Python quantile-regression prediction helper after the R
  `Qtools::rrq()` cleanup; unlike the R path, the Python helper already fits
  `statsmodels` quantile regressions once per requested quantile, so there was
  no duplicate-fit cleanup to retain.
- Re-ran focused current-source cross-platform drift checks after the latest R
  and Python speed checkpoints. RESTAT WaveAll R/Python deterministic outputs
  remained within `5e-5`; RESTAT WaveAll Stata/Python deterministic and
  Gaussian interval matrices remained within `1e-6`, with the largest Gaussian
  interval difference about `4.2e-7`. Multi-illustration remains the noisier
  cross-platform check and continues to show the previously observed fitted and
  interval drift.
- Tested storing Python quantile-regression predictions as numeric `float`
  arrays instead of object arrays. Multi-illustration and RESTAT WaveAll
  snapshots matched at `1e-8`, but single-run timing did not improve, so the
  experiment was not retained.
- Tested batching final Python `scdataMulti()` `Y_pre`/`Y_post` assembly into a
  single concatenation. Common numeric leaves matched at `1e-8`, but the change
  altered the serialized dtype surface and did not improve single-run timing,
  so the experiment was not retained.
- Tested replacing Python `deepcopy()` calls on input pandas DataFrames in
  `scdata()` and `scdataMulti()` with `DataFrame.copy(deep=True)`. Multi-
  illustration and RESTAT WaveAll snapshots matched at `1e-8`, but single-run
  timing did not improve, so the experiment was not retained.
- Tested replacing the balanced-panel `deepcopy(data.unstack().stack(...))` in
  Python `scdata()` with a pandas `.copy()`. Multi-illustration and RESTAT
  WaveAll snapshots matched at `1e-8`, but single-run timing was slower, so the
  experiment was not retained.
- Tested precomputing horizon-level `P @ beta` terms in Python `scpi()` and
  replacing pure-Python sums after ECOS solves with vector dot products.
  Snapshots matched at `1e-8`, but timing was not consistently better across
  multi-illustration and RESTAT WaveAll, so the experiment was not retained.
- Aligned the R and Python multi-unit illustration/benchmark setups with the
  Stata illustration pattern. The updated R/Python multi-illustration
  comparison has no length mismatches; deterministic fitted objects are mostly
  within `5e-5`, with remaining drift concentrated in `b`/`r` and Gaussian
  interval objects. The existing RESTAT WaveAll checks remain the cleaner
  cross-platform validation target.
- Ran the final validation sweep for this modernization round: `R CMD check
  scpi --no-manual --no-build-vignettes` passed with the expected source-prep
  NOTE and restricted-network repository warnings; Python tests passed (`5
  passed`, with expected statsmodels quantile-regression iteration warnings);
  Stata 16 rebuilt/validated `lscpi.mlib`, resolved all six ado commands, and
  regenerated all six Stata help PDFs. The current Stata source tree has no
  `.mata` or `.mo` files, so the rebuilt `lscpi.mlib` contains 0 members.

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
