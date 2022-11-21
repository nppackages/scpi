################################################################################
#' Summary Method for Synthetic Control
#'
#' @description The print method for synthetic control data objects.
#'
#' @param x Class "scdataMulti" object, obtained by calling  \code{\link{scdataMulti}}.
#' @param ... Other arguments.
#'
#' @return No return value, called to print \code{\link{scdataMulti}} results.
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
#' @seealso \code{\link{scdataMulti}} for synthetic control data preparation.
#'
#' Supported methods: \code{\link{print.scdataMulti}}, \code{\link{summary.scdataMulti}}.
#'
#'
#' @export
#'

print.scdataMulti <- function(x, ...) {
    trunits <- length(x$specs$treated.units)
    cat(paste0("Prepared Data for ", trunits, " treated units.\n"))
}

################################################################################
#' Summary Method for Synthetic Control Prediction
#'
#' @description The summary method for synthetic control prediction objects.
#'
#' @param object Class "scdataMulti" object, obtained by calling \code{\link{scdataMulti}}.
#' @param ... Additional arguments
#'
#' @return No return value, called to summarize \code{\link{scdataMulti}} results.
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
#' @seealso \code{\link{scdataMulti}}
#'
#' Supported methods: \code{\link{print.scdataMulti}}, \code{\link{summary.scdataMulti}}.
#'
#' @export

summary.scdataMulti <- function(object, ...) {
    trunits <- object$specs$treated.units

    for (tr in trunits) {
      J       <- object$specs$J[[tr]]
      M       <- object$specs$M[[tr]]
      K       <- object$specs$K[[tr]]
      KM      <- object$specs$KM[[tr]]
      T0      <- object$specs$T0.features[[tr]]

      pt.in   <- object$specs$period.pre[[tr]][1]
      pt.fi   <- object$specs$period.pre[[tr]][length(object$specs$period.pre[[tr]])]

      cat("--------------------------------------------------------------------\n")
      cat(paste0("Synthetic Control - Setup for ", tr, " \n"))
      cat("\n")

      cat(paste("Treated Unit:                              ", tr , "\n", sep = ""))
      cat(paste("Size of the donor pool:                    ", J, "\n", sep = ""))
      cat(paste("Features:                                  ", M, "\n", sep = ""))
      cat(paste("Pre-treatment period:                      ", pt.in, "-", pt.fi, "\n", sep = ""))

      if (M == 1) {
        cat(paste("Pre-treatment periods used in estimation:  ", T0, "\n", sep = ""))
        cat(paste("Covariates used for adjustment:            ", KM, "\n", sep = ""))

      } else {
        cat("Pre-treatment periods used in estimation per feature:\n")
        print(T0)
        cat("Covariates used for adjustment per feature:\n")
        print(K)
      }
      cat("\n")
    }
    cat("--------------------------------------------------------------------\n")

}
