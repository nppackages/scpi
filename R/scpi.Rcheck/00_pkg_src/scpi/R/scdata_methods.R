################################################################################
#' Summary Method for Synthetic Control
#'
#' @description The print method for synthetic control data objects.
#'
#' @param x Class "scdata" object, obtained by calling \code{\link{scdata}}.
#' @param ... Other arguments.
#'
#' @return No return value, called to print \code{\link{scdata}} results.
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
#' @seealso \code{\link{scdata}} for synthetic control data preparation.
#'
#' Supported methods: \code{\link{print.scdata}}, \code{\link{summary.scdata}}.
#'
#' @export
#'
#'

print.scdata <- function(x, ...) {
    tr.unit <- colnames(x$A)
    cat(paste0("Prepared Data for ", tr.unit, ".\n"))
}

################################################################################
#' Summary Method for Synthetic Control Prediction
#'
#' @description The summary method for synthetic control prediction objects.
#'
#' @param object Class "scest" object, obtained by calling \code{\link{scdata}}.
#' @param ... Additional arguments
#'
#' @return No return value, called to summarize \code{\link{scdata}} results.
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
#' @seealso \code{\link{scdata}}
#'
#' Supported methods: \code{\link{print.scdata}}, \code{\link{summary.scdata}}.
#'
#' @export

summary.scdata <- function(object, ...) {

    J       <- object$specs$J
    M       <- object$specs$M
    K       <- object$specs$K
    KM      <- object$specs$KM
    T0      <- object$specs$T0.features
    tr.unit <- colnames(object$A)
    pt.in   <- strsplit(rownames(object$Y.pre)[1], "\\.")[[1]][2]
    pt.fi   <- strsplit(rownames(object$Y.pre)[length(object$Y.pre)], "\\.")[[1]][2]
    pot.in  <- strsplit(rownames(object$Y.post)[1], "\\.")[[1]][2]
    pot.fi  <- strsplit(rownames(object$Y.post)[length(object$Y.post)], "\\.")[[1]][2]
    
    cat("\n")
    cat(paste0("Synthetic Control - Setup\n"))
    cat("\n")

    cat(paste("Treated Unit:                              ", tr.unit, "\n", sep = ""))
    cat(paste("Size of the donor pool:                    ", J, "\n", sep = ""))
    cat(paste("Features:                                  ", M, "\n", sep = ""))
    cat(paste("Pre-treatment period:                      ", pt.in, " || ", pt.fi, "\n", sep = ""))
    cat(paste("Post-treatment period:                     ", pot.in, " || ", pot.fi, "\n", sep = ""))
    
    if (M == 1) {
      cat(paste("Pre-treatment periods used in estimation:  ", T0, "\n", sep = ""))
      cat(paste("Covariates used for adjustment:            ", KM, "\n", sep = ""))

    } else {
      cat("Pre-treatment periods used in estimation per feature:\n")
      print(T0)
      cat("Covariates used for adjustment per feature:\n")
      print(K)
    }

}
