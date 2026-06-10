## Release summary

This is a patch release updating the package to version 4.0.1.

This release:

* updates the CVXR namespace import exclusions for CVXR 1.9.1 compatibility;
* raises the CVXR dependency requirement to `CVXR > 1.9`.

## R CMD check results

Checked locally on Windows 11 x64 using R version 4.6.0 (2026-04-24 ucrt)
with `R CMD check --no-manual`.

0 errors | 0 warnings | 0 notes

Also checked locally with `R CMD check --as-cran` after disabling the CRAN
incoming check because the local R process could not reach the CRAN and
Bioconductor package indexes.

0 errors | 0 warnings | 2 notes

The NOTEs are local-environment artifacts:

* unable to verify the current time remotely;
* a temporary MiKTeX file named `lastMiKTeXException`.

## Reverse dependencies

There are currently no downstream dependencies for this package.
