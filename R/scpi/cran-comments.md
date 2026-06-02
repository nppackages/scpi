## Release summary

This is a major release updating the package to version 4.0.0.

This release:

* updates the maintainer to Matias D. Cattaneo;
* updates package metadata, references, README files, and GitHub workflow
  support;
* adds GitHub issue-tracker metadata through `BugReports`;
* modernizes the Stata workflow for Stata 16 compatibility and distribution
  through a single Mata library;
* adds and validates the Stata `precision(single|double)` option, with
  `double` as the default and `single` preserving the previous behavior;
* improves and validates selected R, Python, and Stata internals while
  preserving numerical behavior.

The previous CRAN maintainer was Filippo Palomba
<fpalomba@princeton.edu>. The current maintainer is Matias D. Cattaneo
<matias.d.cattaneo@gmail.com>. This maintainer change is intentional and has
been coordinated by the package authors.

## R CMD check results

Checked locally on Windows 11 x64 using R version 4.6.0 (2026-04-24 ucrt).

0 errors | 0 warnings | 1 note

The NOTE is the expected CRAN incoming feasibility NOTE for the maintainer
change:

* New maintainer: Matias D. Cattaneo <matias.d.cattaneo@gmail.com>
* Old maintainer: Filippo Palomba <fpalomba@princeton.edu>

## Reverse dependencies

There are currently no downstream dependencies for this package.
