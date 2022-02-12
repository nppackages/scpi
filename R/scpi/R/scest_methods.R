################################################################################
#' Print Method for Synthetic Control Estimation
#'
#' @description The print method for synthetic control estimation objects.
#'
#' @param x Class "scest" object, obtained by calling \code{\link{scest}}.
#' @param ... Other arguments.
#'
#' @author
#' \itemize{
#' \item{Matias Cattaneo, }{Princeton University}
#' \item{Yingjie Feng, }{Tsinghua University}
#' \item{Filippo Palomba, Princeton University (maintainer). \email{fpalomba@princeton.edu}.}
#' \item{Rocio Titiunik, Princeton University}}
#'
#' @seealso \code{\link{scest}} for synthetic control estimation.
#'
#' Supported methods: \code{\link{print.scest}}, \code{\link{summary.scest}}.
#'
#'
#' @export
#' 


print.scest <- function(x, ...) {
  
  Weights    <- round(x$est.results$w, digits = 3)
  if (length(x$est.results$r) > 0) {
    Covariates <- round(x$est.results$r, digits = 3) 
  }
  active.w  <- sum(abs(Weights) > 0)
  
  cat("\n")
  cat("Synthetic Control Estimation - Results\n")
  cat("\n")
  cat(paste("Active donors:", active.w,"\n"))
  cat("\n")
  cat("Coefficients:\n")
  print(cbind(Weights), col.names = F)
  if (length(x$est.results$r) > 0) {
    print(cbind(Covariates), col.names = F)
  }
}




################################################################################
#' Summary Method for Synthetic Control Estimation
#'
#' @description The summary method for synthetic control estimation objects.
#'
#' @param object Class "scest" object, obtained by calling \code{\link{scest}}.
#' @param ... Additional arguments
#'
#' @author
#' \itemize{
#' \item{Matias Cattaneo, }{Princeton University}
#' \item{Yingjie Feng, }{Tsinghua University}
#' \item{Filippo Palomba, Princeton University (maintainer). \email{fpalomba@princeton.edu}.}
#' \item{Rocio Titiunik, Princeton University}}
#' 
#' 
#' @seealso \code{\link{scest}} 
#'
#' Supported methods: \code{\link{print.scest}}, \code{\link{summary.scest}}.
#'
#' @export

summary.scest <- function(object, ...) {
  
  J       <- object$data$specs$J
  M       <- object$data$specs$M
  K       <- object$data$specs$K
  KM      <- object$data$specs$KM
  T0      <- object$data$specs$T0.features
  tr.unit <- colnames(object$data$A)
  pt.in   <- rownames(object$data$Y.pre)[1]
  pt.fi   <- rownames(object$data$Y.pre)[length(object$data$Y.pre)]  
  w.cons  <- object$est.results$w.constr[["name"]]
  if (is.null(object$est.results$w.constr[["Q"]])) {
    w.size <- "-"
  } else {
    w.size  <- round(object$est.results$w.constr[["Q"]], 3)
  }
  cat("\n")
  cat("Synthetic Control Estimation - Setup\n")
  cat("\n")
  
  cat(paste("Constraint Type:                           ", w.cons, "\n", sep = ""))
  cat(paste("Constraint Size (Q):                       ", w.size, "\n", sep = ""))
  cat(paste("Treated Unit:                              ", tr.unit,"\n", sep = ""))
  cat(paste("Size of the donor pool:                    ", J,"\n", sep = ""))
  cat(paste("Features:                                  ", M,"\n", sep = ""))
  cat(paste("Pre-treatment period:                      ", pt.in,"-",pt.fi,"\n", sep = ""))
  
  if (M == 1) {
    cat(paste("Pre-treatment periods used in estimation:  ",T0,"\n", sep = ""))
    cat(paste("Covariates used for adjustment:            ",KM,"\n", sep = ""))
    
  } else {
    cat("Pre-treatment periods used in estimation per feature:\n")
    print(T0)
    cat("Covariates used for adjustment per feature:\n")
    print(K)
  }
  
  print(object)
}