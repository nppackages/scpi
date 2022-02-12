#' @title scpi: A package to compute Synthetic Control Prediction Intervals
#'
#' @description The package computes point estimates and prediction intervals for Synthetic Control methods.
#'
#' @author
#' \itemize{
#' \item{Matias Cattaneo, }{Princeton University}
#' \item{Yingjie Feng, }{Tsinghua University}
#' \item{Filippo Palomba, Princeton University (maintainer). \email{fpalomba@princeton.edu}.}
#' \item{Rocio Titiunik, Princeton University}}
#' 
#' @references
#' \itemize{
#' \item{\href{https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, M. D., Feng, Y., & Titiunik, R. 
#' (2021)}. Prediction intervals for synthetic control methods. \emph{Journal of the American Statistical Association}, 116(536), 1865-1880.}
#' \item{\href{https://nppackages.github.io/references/Cattaneo-Feng-Palomba-Titiunik_2022_scpi.pdf}{Cattaneo, M. D., Feng, Y., Palomba F., and Titiunik, R. (2022)}}
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
#' 
#' @aliases scpi-package
"_PACKAGE"

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
