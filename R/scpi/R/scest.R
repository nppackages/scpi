###############################################################################

#' @title Estimation of Synthetic Control
#'
#' @description The command implements estimation procedures for Synthetic Control (SC) methods using least square, lasso, ridge, or simplex-type constraints according to
#'  \href{https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, M. D., Feng, Y., & Titiunik, R. (2021)}. 
#'
#' Companion \href{https://www.stata.com/}{Stata} and \href{https://www.python.org/}{Python} packages are described in \href{https://arxiv.org/abs/2202.05984}{Cattaneo, Feng, Palomba, and Titiunik (2022)}.
#'
#' Companion commands are: \link{scdata} for data preparation, \link{scpi} for inference procedures, and \link{scplot} for plots.
#' 
#' Related Stata, R, and Python packages useful for inference in SC designs are described in the following website:
#' 
#' \href{ https://nppackages.github.io/scpi/}{ https://nppackages.github.io/scpi/}
#' 
#' For an introduction to synthetic control methods, see \href{https://economics.mit.edu/files/17847}{Abadie (2021)} and references therein.
#' 
#' @param data a class `scpi_data' object, obtained by calling \code{\link{scdata}}.
#' @param w.constr a list specifying the constraint set the estimated weights of the donors must belong to.
#' \code{w.constr} can contain up to four objects:
#' - `\code{p}', a string indicating the norm to be constrained (\code{p} should be one of "no norm", "L1", and "L2")
#' - `\code{dir}', a string indicating whether the constraint on the norm is an equality ("==") or inequality ("<=")
#' - `\code{Q}', a scalar defining the value of the constraint on the norm
#' - `\code{lb}', a scalar defining the lower bound on the weights. It can be either 0 or \code{-Inf}.
#' - `\code{name}', a character selecting one of the default proposals.
#' See the \strong{Details} section for more.
#' @param V a weighting matrix to be used when minimizing the sum of squared
#' residuals
#' \deqn{(\mathbf{A}-\mathbf{B}\mathbf{w}-\mathbf{C}\mathbf{r})'\mathbf{V}(\mathbf{A}-\mathbf{B}\mathbf{w}-\mathbf{C}\mathbf{r})}
#' The default is the identity matrix, so equal weight is given to all observations.
#'
#' @param opt.list a list specifying the stopping criteria and the algorithm for the underlying optimizer (\code{\link{nloptr}} or \code{\link{CVXR}}) for point estimation.
#' See the \strong{Details} section for more.
#' 
#' @param plot a logical specifying whether \code{\link{scplot}} should be called and a plot saved in the current working directory. For more options see \code{\link{scplot}}.
#' @param plot.name a string containing the name of the plot (the format is by default .png). For more options see \code{\link{scplot}}.
#' @param plot.path a string containing the path at which the plot should be saved (default is output of \code{getwd()}.)
#' @param save.data a character specifying the name and the path of the saved dataframe containing the processed data used to produce the plot. 
#' 
#' @return
#' The function returns an object of class 'scpi_scest' containing two lists. The first list is labeled 'data' and contains used data as returned by \code{\link{scdata}} and some other values.
#' \item{A}{a matrix containing pre-treatment features of the treated unit.}
#' \item{B}{a matrix containing pre-treatment features of the control units.}
#' \item{C}{a matrix containing covariates for adjustment.}
#' \item{P}{a matrix whose rows are the vectors used to predict the out-of-sample series for the synthetic unit.}
#' \item{Y.pre}{a matrix containing the pre-treatment outcome of the treated unit.}
#' \item{Y.post}{a matrix containing the post-treatment outcome of the treated unit.}
#' \item{Y.donors}{a matrix containing the pre-treatment outcome of the control units.}
#' \item{specs}{a list containing some specifics of the data:
#' \itemize{
#' \item{\code{J}, the number of control units}
#' \item{\code{K}, a numeric vector with the number of covariates used for adjustment for each feature}
#' \item{\code{KM}, the total number of covariates used for adjustment}
#' \item{\code{M}, number of features}
#' \item{\code{period.pre}, a numeric vector with the pre-treatment period}
#' \item{\code{period.post}, a numeric vector with the post-treatment period}
#' \item{\code{T0.features}, a numeric vector with the number of periods used in estimation for each feature}
#' \item{\code{T1.outcome}, the number of post-treatment periods}
#' \item{\code{glob.cons}, for internal use only}
#' \item{\code{out.in.features}, for internal use only}}}
#' 
#' The second list is labeled 'est.results' and contains estimation results.
#' \item{w}{a matrix containing the estimated weights of the donors.}
#' \item{r}{a matrix containing the values of the covariates used for adjustment.}
#' \item{b}{a matrix containing \eqn{\mathbf{w}} and \eqn{\mathbf{r}}.}
#' \item{Y.pre.fit}{a matrix containing the estimated pre-treatment outcome of the SC unit.}
#' \item{Y.post.fit}{a matrix containing the estimated post-treatment outcome of the SC unit.}
#' \item{A.hat}{a matrix containing the predicted values of the features of the treated unit.}
#' \item{res}{a matrix containing the residuals \eqn{\mathbf{A}-\widehat{\mathbf{A}}}.}
#' \item{V}{a matrix containing the weighting matrix used in estimation.}
#' \item{w.constr}{a list containing the specifics of the constraint set used on the weights.}
#' 
#' 
#' @details
#' \itemize{
#' \item{\strong{Estimation of Weights.} \code{w.constr} specifies the constraint set on the weights. First, the element
#' \code{p} allows the user to choose between imposing a constraint on either the L1 (\code{p = "L1"})
#' or the L2 (\code{p = "L2"}) norm of the weights and imposing no constraint on the norm (\code{p = "no norm"}).
#' Second, \code{Q} specifies the value of the constraint on the norm of the weights.
#' Third, \code{lb} sets the lower bound of each component of the vector of weights.
#' Fourth, \code{dir} sets the direction of the constraint on the norm in case \code{p = "L1"}
#' or \code{p = "L2"}. If \code{dir = "=="}, then
#' \deqn{||\mathbf{w}||_p = Q,\:\:\: w_j \geq lb,\:\: j =1,\ldots,J}
#' If instead \code{dir = "<="}, then
#' \deqn{||\mathbf{w}||_p \leq Q,\:\:\: w_j \geq lb,\:\: j =1,\ldots,J}
#' If instead \code{dir = "NULL"} no constraint on the norm of the weights is imposed.
#'
#' An alternative to specifying an ad-hoc constraint set on the weights would be
#' choosing among some popular types of constraints. This can be done by including the element
#' `\code{name}' in the list \code{w.constr}. The following are available options:
#' \itemize{
#' \item {If \code{name == "simplex"} (the default), then
#' \deqn{||\mathbf{w}||_1 = 1,\:\:\: w_j \geq 0,\:\: j =1,\ldots,J.}}
#' 
#' \item {If \code{name == "lasso"}, then
#' \deqn{||\mathbf{w}||_1 \leq Q,}
#' where \code{Q} is by default equal to 1 but it can be provided as an element of the list (eg. \code{w.constr =
#' list(name = "lasso", Q = 2)}).}
#' 
#' \item{If \code{name == "ridge"}, then
#' \deqn{||\mathbf{w}||_2 \leq Q,}
#' where \code{Q} is a tuning parameter that is by default computed as 
#' \deqn{(J+KM) \widehat{\sigma}_u^{2}/||\widehat{\mathbf{w}}_{OLS}||_{2}^{2}}
#' where \eqn{J} is the number of donors and \eqn{KM} is the total number of covariates used for adjustment.
#' The user can provide \code{Q} as an element of the list (eg. \code{w.constr =
#' list(name = "ridge", Q = 1)}).}
#'
#' \item{If \code{name == "ols"}, then the problem is unconstrained and the vector of weights
#' is estimated via ordinary least squares.}
#' }
#' }
#' 
#' \item{\strong{Algorithm Options.} The default is a sequential quadratic programming (SQP) algorithm for nonlinearly constrained gradient-based optimization 
#' (\code{algorithm = 'NLOPTR_LD_SLSQP'}) for all cases not involving the L1 norm. 
#' For a complete list of algorithms see \href{https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/}{the official \code{nloptr} documentation}.
#' The other default values are \code{maxeval = 5000}, \code{ftol_res = 1.0e-8}, \code{ftol_abs = 1.0e-8}, \code{xtol_res = 1.0e-8}, \code{xtol_abs = 1.0e-8},
#' \code{tol_constraints_eq = 1.0e-8}, and, finally, \code{tol_constraints_ineq = 1.0e-8}. More information on the stopping criteria can be obtained running
#' \code{nloptr.print.options()} or \code{nloptr.get.default.options()}. f the optimization involves the L1 norm then \code{CVXR} is used for optimization.  
#' More information on the stopping criteria can be obtained reading \href{https://cvxr.rbind.io/}{the official documentation}. }
#' }
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
#' @seealso \code{\link{scdata}}, \code{\link{scpi}}, \code{\link{scplot}}
#'
#' @examples
#' 
#' data <- scpi_germany
#' 
#' df <- scdata(df = data, id.var = "country", time.var = "year", 
#'              outcome.var = "gdp", period.pre = (1960:1990), 
#'              period.post = (1991:2003), unit.tr = "West Germany",
#'              unit.co = unique(data$country)[-7], constant = TRUE,
#'              cointegrated.data = TRUE)
#'              
#' result <- scest(df, w.constr = list(name = "simplex", Q = 1))
#' result <- scest(df, w.constr = list(lb = 0, dir = "==", p = "L1", Q = 1))
#' 
#' 
#' @export

scest <- function(data,
                  w.constr  = NULL,
                  V         = NULL,
                  opt.list  = NULL,
                  plot      = FALSE,
                  plot.name = NULL,
                  plot.path = NULL,
                  save.data = NULL) {

  ##########################
  if (class(data)[1] != "scpi_data") {
    stop("data should be the object returned by running scdata!")
  }

  if (is.null(w.constr) == FALSE) { # The user has specified W

    if (is.list(w.constr) == FALSE) {
      stop("w.constr should be a list!")
    }

    if (!"name" %in% names(w.constr)) {
      w.constr[["name"]] <- "NULL"
    } else {
      if (!w.constr[["name"]] %in% c("simplex","lasso","ridge","ols")) {
        stop("If 'name' is specified in w.constr, then it should be chosen among
             'simplex', 'lasso', 'ridge', 'ols'.")
      }
    }
  }
  
  if (is.null(opt.list) == FALSE) {
    if (is.list(opt.list) == FALSE) {
      stop("The object opt.list should be a list!")
    }
  }  
  
  # Data matrices
  A <- data$A
  B <- data$B
  C <- data$C
  P <- data$P
  Z <- cbind(B, C)
  Y.donors <- data$Y.donors

  # Data specs
  J               <- data$specs$J
  KM              <- data$specs$KM
  K               <- data$specs$K
  M               <- data$specs$M
  T0.features     <- data$specs$T0.features
  out.in.features <- data$specs$out.in.features

  ##########################
  ## Set up the estimation problem

  # Create weighting matrix
  if (is.null(V) == FALSE){
    if (is.matrix(V) == FALSE) {
      stop("The object V should be a matrix!")
    }
    
    if (nrow(V) != nrow(B) | ncol(V) != nrow(B)) {
      stop(paste0("V should be a ", nrow(B), "x", nrow(B)," matrix, but currently it is
                  a ", nrow(V),"x", ncol(V)," matrix!"))
    }
  }
  
  if (is.null(V)) {
    V <- diag(dim(B)[1])
  }

  
  # Create constraints
  w.constr <- w.constr.OBJ(w.constr,A,Z,V,J,KM)
  
  # if (!length(w.constr$lb) %in% c(1,J)) {
  #   stop("The length of 'lb' should be either one or the number of donors.")
  # }

  # Estimate SC
  b <- b.est(A = A, Z = Z, J = J, KM = KM, w.constr = w.constr, V = V, QQ = w.constr[["Q"]], opt.list = opt.list)

  ##########################
  ## Create useful objects  
  # Store point estimates
  if (KM == 0) {
    w <- b                       # Estimated weights
    r <- NULL                    # Loading factors of additional covariates
  } else {
    w <- b[1:J]
    r <- b[(J+1):(J+KM)]
  }

  # Fitted values and residuals
  A.hat    <- Z %*% b                            # Pre-treatment fit of synthetic unit
  res      <- A  -  A.hat                        # Estimation residual u
  
  # Pre-treatment fit of outcome of interest
  if (out.in.features == TRUE) {
    fit.pre  <- A.hat[1:T0.features[1], , drop= FALSE]
    names    <- strsplit(rownames(fit.pre), "\\.")
    rownames(fit.pre) <- unlist(lapply(names, "[[", 2))

  } else {
    fit.pre  <- Y.donors %*% w
  }
  # Post-treatment prediction of outcome of interest
  fit.post <- P %*% b             

  est.results <- list(b = b,
                      w = w,
                      r = r,
                      Y.pre.fit = fit.pre,
                      Y.post.fit = fit.post,
                      A.hat = A.hat,
                      res = res,
                      V = V,
                      w.constr = w.constr)


  df   <- list(A = data$A,
               B = data$B,
               C = data$C,
               P = data$P,
               Z = Z,
               specs  = data$specs,
               Y.pre  = data$Y.pre,
               Y.post = data$Y.post)

  ##########################
  ## Return to the user
  to.return <- list(data = df, est.results = est.results)
  class(to.return) <- 'scest'

  ##################################################
  ## Plot
  if (plot == TRUE) {
    if (is.null(plot.name) == FALSE) {
      fig.name <- plot.name
    } else {
      fig.name <- "scest_default_plot"
    }
    
    if (is.null(plot.path) == FALSE) {
      fig.path <- plot.path
    } else {
      fig.path <- getwd()
    }
    
    scplot(result = to.return, fig.path = fig.path,
           fig.name = fig.name, fig.format = "png", save.data = save.data)

  }



  return(to.return)
}

###############################################################################
