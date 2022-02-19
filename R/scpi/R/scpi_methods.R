################################################################################
#' Print Method for Synthetic Control Inference
#'
#' @description The print method for synthetic control inference objects.
#'
#' @param x Class "scpi" object, obtained by calling \code{\link{scpi}}.
#' @param ... Other arguments.
#'
#' @return No return value, called to print \code{\link{scpi}} results.
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
#' @seealso \code{\link{scpi}} for synthetic control inference
#'
#' Supported methods: \code{\link{print.scpi}}, \code{\link{summary.scpi}}.
#'
#'
#' @export
#' 


print.scpi <- function(x, ...) {
  
  args <- list(...)
  if (is.null(args[['already.est']])) {print_est <- TRUE} else {print_est <- FALSE}
  
  e.method  <- x$inference.results$e.method
  Y.tr.post <- x$data$Y.post
  Y.sc.post <- round(x$est.results$Y.post.fit, 3)
  colnames(Y.tr.post) <- "Treated"
  colnames(Y.sc.post) <- "Synthetic"
  
  if (e.method == "gaussian") {
    CI <- x$inference.results$CI.all.gaussian[, 1:2, drop = FALSE]
    
  } else if (e.method == "ls") {
    CI <- x$inference.results$CI.all.ls[, 1:2, drop = FALSE]
    
  } else if (e.method == "qreg") {
    CI <- x$inference.results$CI.all.qreg[, 1:2, drop = FALSE]
    
  } else if (e.method == "all") {
    CI <- cbind(x$inference.results$CI.all.gaussian[, 1:2, drop = FALSE], 
                x$inference.results$CI.all.ls[, 1:2, drop = FALSE], 
                x$inference.results$CI.all.qreg[, 1:2, drop = FALSE])
  }
  
  
  if (print_est == TRUE) {  # If false avoids printing twice when summmary is called
    xx        <- x
    class(xx) <- 'scest'
    print(xx)
  }
  
  
  cat("\n")
  cat("Synthetic Control Inference - Results\n")
  cat("\n")
  if (e.method == "gaussian") {
    cat("  Inference with subgaussian bounds   \n")
    
  } else if (e.method == "ls") {
    cat("  Inference with location-scale model   \n")
    
  } else if (e.method == "qreg") {
    cat("  Inference with quantile regression  \n")
    
  } else if (e.method == "all") {
    cat("                             Subgaussian            Location Scale        Quantile Reg   \n")  
  }
  
  # Check for eventual missing values in the donor pool that didn't allow
  # inference to be conducted for some periods
  inf.con <- rownames(Y.sc.post)
  print(cbind(Y.tr.post[inf.con,,drop = FALSE], Y.sc.post, round(CI, digits = 3)), ...)
}





################################################################################
#' Summary Method for Synthetic Control Inference
#'
#' @description The summary method for synthetic control inference objects.
#'
#' @param object Class "scpi" object, obtained by calling \code{\link{scpi}}.
#' @param ... Additional arguments
#'
#' @return No return value, called to summarize \code{\link{scpi}} results.
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
#' 
#' @seealso \code{\link{scpi}} 
#'
#' Supported methods: \code{\link{print.scpi}}, \code{\link{summary.scpi}}.
#'
#' @export


summary.scpi <- function(object, ...) {
  xx        <- object
  class(xx) <- 'scest'
  summary(xx)
  
  
  cat("\n")
  cat("Synthetic Control Inference - Setup\n")
  cat("\n")
  
  if (object$inference.results$u.user == FALSE) {
    cat(paste("In-sample Inference:                       ", "\n", sep = ""))
    cat(paste("     Misspecified model                    ", object$inference.results$u.missp, "\n", sep = ""))
    cat(paste("     Order of polynomial (B)               ", object$inference.results$u.order,"\n", sep = ""))
    cat(paste("     Lags (B)                              ", object$inference.results$u.lags,  "\n", sep = ""))
    cat(paste("     Variance-Covariance Estimator         ", object$inference.results$u.sigma, "\n", sep = ""))  
  } else {
    cat(paste("In-sample Inference:                       ", "\n", sep = ""))
    cat("      User provided \n")    
  }
  
  cat("\n")
  if (object$inference.results$e.user == FALSE) {
    cat(paste("Out-of-sample Inference:                   ", "\n", sep = ""))
    cat(paste("     Method                                ", object$inference.results$e.method, "\n", sep = ""))
    cat(paste("     Order of polynomial (B)               ", object$inference.results$u.order,  "\n", sep = ""))
    cat(paste("     Lags (B)                              ", object$inference.results$u.lags,   "\n", sep = ""))
  } else {
    cat(paste("Out-of-sample Inference:                   ", "\n", sep = ""))
    cat("      User provided \n")    
  }
  
  cat("\n")
  
  print(object, already.est = TRUE, ...)
  
  if (object$est.results$w.constr[["name"]] == "ridge") {
    warning("The current version of the package does not take into account a small probability loss in prediction intervals when using 
            'ridge' as estimation method! An updated version will be rolled out soon!!")
  }
}