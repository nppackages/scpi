#' @title scpi: A package to compute Synthetic Control Prediction Intervals
#'
#' @description The package implements estimation, inference procedures, and produce plots for Synthetic Control (SC) methods using least square, lasso, ridge, or simplex-type 
#' constraints according to \href{https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, M. D., Feng, Y., & Titiunik, R. (2021)}.
#' 
#' Included functions are: \link{scdata} for data preparation, \link{scest} for point estimation, \link{scpi} for inference procedures, and \link{scplot} for plots.
#' 
#' \code{print()} and \code{summary()} methods are available for \code{\link{scest}} and \code{\link{scpi}}.
#'
#' Companion \href{https://www.stata.com/}{Stata} and \href{https://www.python.org/}{Python} packages are described in \href{https://arxiv.org/abs/2202.05984}{Cattaneo, Feng, Palomba, and Titiunik (2022)}.
#' 
#' Related Stata, R, and Python packages useful for inference in SC designs are described in the following website:
#' 
#' \href{ https://nppackages.github.io/scpi/}{ https://nppackages.github.io/scpi/}
#' 
#' For an introduction to synthetic control methods, see \href{https://economics.mit.edu/files/17847}{Abadie (2021)} and references therein.
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
#' \itemize{
#' \item{\href{https://economics.mit.edu/files/17847}{Abadie, A. (2021)}. Using synthetic controls: Feasibility, data requirements, and methodological aspects.
#' \emph{Journal of Economic Literature}, 59(2), 391-425.}
#' \item{\href{https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, M. D., Feng, Y., & Titiunik, R. 
#' (2021)}. Prediction intervals for synthetic control methods. \emph{Journal of the American Statistical Association}, 116(536), 1865-1880.}
#' \item{\href{https://arxiv.org/abs/2202.05984}{Cattaneo, M. D., Feng, Y., Palomba F., and Titiunik, R. (2022).}
#' scpi: Uncertainty Quantification for Synthetic Control Estimators, \emph{arXiv}:2202.05984.}
#' }
#'
#' 
#' @importFrom dplyr lag
#' @importFrom dplyr left_join
#' @importFrom ECOSolveR ECOS_csolve
#' @importFrom fastDummies dummy_cols
#' @importFrom magrittr %>%
#' @importFrom Matrix bdiag
#' @importFrom purrr map
#' @importFrom Qtools rrq
#' @importFrom stringr str_remove
#' @importFrom stringr str_split
#' @importFrom tibble is_tibble
#' @importFrom utils flush.console
#' 
#' @import abind
#' @import CVXR
#' @import doSNOW
#' @import doRNG
#' @import foreach
#' @import ggplot2
#' @import nloptr
#' @import parallel
#' @import tidyr
#' 
#' @rawNamespace import(stats, except = c(lag,filter,power))
#' @rawNamespace import(rlang, except = c(is_vector,is_complex))
#' 
#' @aliases scpi-package
"_PACKAGE"

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
