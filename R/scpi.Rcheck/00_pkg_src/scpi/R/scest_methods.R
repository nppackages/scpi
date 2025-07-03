################################################################################
#' Print Method for Synthetic Control Methods
#'
#' @description The print method for synthetic control prediction fitted objects.
#'
#' @param x Class "scest" object, obtained by calling \code{\link{scest}}.
#' @param ... Other arguments.
#'
#' @return No return value, called to print \code{\link{scest}} results.
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
#' @seealso \code{\link{scest}} for synthetic control prediction.
#'
#' Supported methods: \code{\link{print.scest}}, \code{\link{summary.scest}}, \code{\link{coef.scest}}.
#'
#'
#' @export
#'


print.scest <- function(x, ...) {
  I <- x$data$specs$I

  if (I == 1) {
    Weights    <- round(x$est.results$w, digits = 3)

    if (length(x$est.results$r) > 0) {
      Covariates <- round(x$est.results$r, digits = 3)
    }
    active.w  <- sum(abs(Weights) > 0)

    names <- strsplit(names(Weights), "\\.")
    names <- unlist(lapply(names, "[[", 2))  # Get all control units
    names(Weights) <- names

    cat("\n")
    cat("Synthetic Control Prediction - Results\n")
    cat("\n")
    cat(paste("Active donors:", active.w,"\n"))
    cat("\n")
    cat("Coefficients:\n")
    print(cbind(Weights), col.names = FALSE)
    if (length(x$est.results$r) > 0) {
      print(cbind(Covariates), col.names = FALSE)
    }
  } else if (I > 1) {
    treated.units <- x$data$specs$treated.units
    Weights <- round(x$est.results$w, digits = 3)
    W.list <- mat2list(as.matrix(Weights))

    names <- strsplit(names(Weights), "\\.")
    names <- unlist(lapply(names, "[[", 2))  # Get all control units
    co.units <- unique(names)
    to.print <- data.frame(control.unit = co.units)

    for (i in seq_len(I)) {
      names <- strsplit(rownames(W.list[[i]]), "\\.")
      names <- unlist(lapply(names, "[[", 2))  # Get all control units
      to.merge <- data.frame(control.unit = names, W.list[[i]])
      names(to.merge) <- c("control.unit", treated.units[i])
      to.print <- merge(to.print, to.merge, by="control.unit", all = TRUE)
    }

    mat.print <- as.matrix(to.print[,-1, drop = FALSE])
    rownames(mat.print) <- to.print[,1, drop = TRUE]

    active.w <- colSums(abs(mat.print) > 0, na.rm = TRUE)

    cat("\n")
    cat("Synthetic Control Prediction - Results\n")
    cat("\n")
    cat(paste("Active donors:\n"))
    print(active.w)
    cat("\n")
    cat("Coefficients:\n")
    print(mat.print)
  }
}




################################################################################
#' Summary Method for Synthetic Control Prediction
#'
#' @description The summary method for synthetic control prediction fitted objects.
#'
#' @param object Class "scest" object, obtained by calling \code{\link{scest}}.
#' @param ... Additional arguments
#'
#' @return No return value, called to summarize \code{\link{scest}} results.
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
#' @seealso \code{\link{scest}}
#'
#' Supported methods: \code{\link{print.scest}}, \code{\link{summary.scest}}, \code{\link{coef.scest}}.
#'
#' @export

summary.scest <- function(object, ...) {

  if (object$data$specs$I == 1) {
    J       <- object$data$specs$J
    M       <- object$data$specs$M
    K       <- object$data$specs$K
    KM      <- object$data$specs$KM
    T0      <- object$data$specs$T0.features
    tr.unit <- colnames(object$data$A)
    pt.in   <- strsplit(rownames(object$data$Y.pre)[1], "\\.")[[1]][2]
    pt.fi   <- strsplit(rownames(object$data$Y.pre)[length(object$data$Y.pre)], "\\.")[[1]][2]  
    w.cons  <- object$est.results$w.constr[["name"]]
    if (is.null(object$est.results$w.constr[["Q"]])) {
      w.size <- "-"
    } else {
      w.size  <- round(object$est.results$w.constr[["Q"]], 3)
    }
    cat("\n")
    cat(paste0("Synthetic Control Prediction - Setup\n"))
    cat("\n")

    cat(paste("Constraint Type:                           ", w.cons, "\n", sep = ""))
    cat(paste("Constraint Size (Q):                       ", w.size, "\n", sep = ""))
    cat(paste("Treated Unit:                              ", tr.unit,"\n", sep = ""))
    cat(paste("Size of the donor pool:                    ", J,"\n", sep = ""))
    cat(paste("Features:                                  ", M,"\n", sep = ""))
    cat(paste("Pre-treatment period:                      ", pt.in,"-",pt.fi,"\n", sep = ""))

    if (M == 1) {
      cat(paste("Pre-treatment periods used in prediction:  ",T0,"\n", sep = ""))
      cat(paste("Covariates used for adjustment:            ",KM,"\n", sep = ""))

    } else {
      cat("Pre-treatment periods used in prediction per feature:\n")
      print(T0)
      cat("Covariates used for adjustment per feature:\n")
      print(K)
    }

  } else {
    trunits <- object$data$specs$treated.units

    for (tr in trunits) {
      J       <- object$data$specs$J[[tr]]
      M       <- object$data$specs$M[[tr]]
      K       <- object$data$specs$K[[tr]]
      KM      <- object$data$specs$KM[[tr]]
      T0      <- object$data$specs$T0.features[[tr]]

      pt.in   <- object$data$specs$period.pre[[tr]][1]
      pt.fi   <- object$data$specs$period.pre[[tr]][length(object$data$specs$period.pre[[tr]])]
      w.cons  <- object$est.results$w.constr[[tr]][["name"]]
      if (is.null(object$est.results$w.constr[[tr]][["Q"]])) {
        w.size <- "-"
      } else {
        w.size  <- round(object$est.results$w.constr[[tr]][["Q"]], 3)
      }
      cat("--------------------------------------------------------------------\n")
      cat(paste0("Synthetic Control Prediction - Setup for ",tr," \n"))
      cat("\n")

      cat(paste("Constraint Type:                           ", w.cons, "\n", sep = ""))
      cat(paste("Constraint Size (Q):                       ", w.size, "\n", sep = ""))
      cat(paste("Treated Unit:                              ", tr , "\n", sep = ""))
      cat(paste("Size of the donor pool:                    ", J, "\n", sep = ""))
      cat(paste("Features:                                  ", M, "\n", sep = ""))
      cat(paste("Pre-treatment period:                      ", pt.in, "-", pt.fi, "\n", sep = ""))

      if (M == 1) {
        cat(paste("Pre-treatment periods used in prediction:  ", T0, "\n", sep = ""))
        cat(paste("Covariates used for adjustment:            ", KM, "\n", sep = ""))

      } else {
        cat("Pre-treatment periods used in prediction per feature:\n")
        print(T0)
        cat("Covariates used for adjustment per feature:\n")
        print(K)
      }
      cat("\n")
    }
    cat("--------------------------------------------------------------------\n")
 }

  print.scest(object)
}


################################################################################
#' Coef Method for Synthetic Control Methods
#'
#' @description The coef method for synthetic control prediction fitted objects.
#'
#' @param object Class "scest" object, obtained by calling \code{\link{scest}}.
#' @param ... Other arguments (eg. \code{ncols}).
#'
#' @return No return value, called to show \code{\link{scest}} constructed weights.
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
#' @seealso \code{\link{scest}} for synthetic control prediction.
#'
#' Supported methods: \code{\link{print.scest}}, \code{\link{summary.scest}}, \code{\link{coef.scest}}.
#'
#'
#' @export
#'

coef.scest <- function(object, ...) {

  args <- list(...)
  I <- object$data$specs$I

  if (is.null(args[['ncols']])) {
    ncols <- 1 * (I <= 3) + 1 * (I > 3) + 1 * (I > 6)
  } else {
    ncols <- args[['ncols']]
  }

  w <- object$est.results$w
  aux <- data.frame(w = w,
                    donor = unlist(lapply(strsplit(names(w), "\\."), "[[", 2)),
                    treated = unlist(lapply(strsplit(names(w), "\\."), "[[", 1)))

  ggplot() +
    geom_point(data = aux, aes(x = .data$donor, y = .data$w, size = abs(.data$w))) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    facet_wrap(~treated, ncol = ncols) +
    xlab("") + ylab("Weight") +
    theme(legend.position = "none",
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2, size = 10),
          axis.text.y = element_text(size=13),
          axis.title.x = element_text(size = 13),
          axis.title.y = element_text(size = 13),
          panel.grid.minor = element_blank())

}
