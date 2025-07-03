#' @title \code{scpi}: A Package to Compute Synthetic Control Prediction Intervals With Multiple Treated Units and Staggered Adoption
#'
#' @description The package implements estimation, inference procedures, and produces plots for Synthetic Control (SC) methods
#' using least squares, lasso, ridge, or simplex-type
#' constraints. Uncertainty is quantified using prediction intervals according to
#' \insertCite{cattaneo2021methodological-JASA;textual}{scpi} and \insertCite{cattaneo2025methodological-RESTAT;textual}{scpi}.
#'
#' Included functions are: \link{scdata} and \link{scdataMulti} for data preparation, \link{scest} for point estimation,
#' \link{scpi} for inference procedures, and \link{scplot} and \link{scplotMulti} for plots.
#'
#' \code{print()} and \code{summary()} methods are available for \code{\link{scest}} and \code{\link{scpi}}.
#'
#' Companion \href{https://www.stata.com/}{Stata} and \href{https://www.python.org/}{Python} packages are described in
#' \insertCite{cattaneo2025software-JSS;textual}{scpi}.
#'
#' Related Stata, R, and Python packages useful for inference in SC designs are described in the following website:
#'
#' \href{ https://nppackages.github.io/scpi/}{ https://nppackages.github.io/scpi/}
#'
#' For an introduction to synthetic control methods, see \insertCite{abadie2021UsingSyntheticControls;textual}{scpi} and references therein.
#'
#' @author
#' Matias Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.
#'
#' Yingjie Feng, Tsinghua University. \email{fengyj@sem.tsinghua.edu.cn}.
#'
#' Filippo Palomba, Princeton University (maintainer). \email{fpalomba@princeton.edu}.
#'
#' Rocio Titiunik, Princeton University. \email{titiunik@princeton.edu}.
#'
#' @references
#'  \insertAllCited{}
#'
#' @importFrom dplyr lag
#' @importFrom dplyr left_join
#' @importFrom fastDummies dummy_cols
#' @importFrom magrittr %>%
#' @importFrom Matrix bdiag
#' @importFrom MASS ginv
#' @importFrom methods as
#' @importFrom methods is
#' @importFrom purrr map
#' @importFrom Rdpack reprompt
#' @importFrom reshape2 melt
#' @importFrom Qtools rrq
#' @importFrom stringr str_remove
#' @importFrom stringr str_split
#' @importFrom tibble is_tibble
#' @importFrom utils flush.console
#'
#' @import abind
#' @import CVXR
#' @import doSNOW
#' @import ECOSolveR
#' @import foreach
#' @import ggplot2
#' @import parallel
#' @import tidyr
#'
#' @rawNamespace import(stats, except = c(lag,filter,power))
#' @rawNamespace import(rlang, except = c(is_vector,is_complex))
#'
#' @aliases scpi-package
"_PACKAGE"

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1")  utils::globalVariables(c("."))
