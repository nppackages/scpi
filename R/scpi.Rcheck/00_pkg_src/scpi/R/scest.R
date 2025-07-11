###############################################################################

#' @title Prediction of Synthetic Control
#'
#' @description The command implements estimation procedures for Synthetic Control (SC) methods using least squares, lasso,
#' ridge, or simplex-type constraints. For more information see
#' \insertCite{cattaneo2021methodological-JASA;textual}{scpi} and
#' \insertCite{cattaneo2025methodological-RESTAT;textual}{scpi}.
#'
#' Companion \href{https://www.stata.com/}{Stata} and \href{https://www.python.org/}{Python} packages are described in
#' \insertCite{cattaneo2025software-JSS;textual}{scpi}.
#'
#' Companion commands are: \link{scdata} and \link{scdataMulti} for data preparation in the single and multiple treated unit(s) cases, respectively,
#' \link{scpi} for inference procedures, \link{scplot} and \link{scplotMulti} for plots in the single and multiple treated unit(s) cases, respectively.
#'
#' Related Stata, R, and Python packages useful for inference in SC designs are described in the following website:
#'
#' \href{ https://nppackages.github.io/scpi/}{ https://nppackages.github.io/scpi/}
#'
#' For an introduction to synthetic control methods, see \insertCite{abadie2021UsingSyntheticControls;textual}{scpi} and references therein.
#'
#' @param data a class 'scdata' object, obtained by calling \code{\link{scdata}}, or class 'scdataMulti' obtained via \code{\link{scdataMulti}}.
#' @param w.constr a list specifying the constraint set the estimated weights of the donors must belong to.
#' \code{w.constr} can contain up to four objects:
#' - `\code{p}', a string indicating the norm to be constrained (\code{p} should be one of "no norm", "L1", and "L2")
#' - `\code{dir}', a string indicating whether the constraint on the norm is an equality ("==") or inequality ("<=")
#' - `\code{Q}', a scalar defining the value of the constraint on the norm
#' - `\code{lb}', a scalar defining the lower bound on the weights. It can be either 0 or \code{-Inf}.
#' - `\code{name}', a character selecting one of the default proposals.
#' See the \strong{Details} section for more.
#' @param V specifies the type of weighting matrix to be used when minimizing the sum of squared residuals
#' \deqn{(\mathbf{A}-\mathbf{B}\mathbf{w}-\mathbf{C}\mathbf{r})'\mathbf{V}(\mathbf{A}-\mathbf{B}\mathbf{w}-\mathbf{C}\mathbf{r})}
#' The default is the identity matrix, so equal weight is given to all observations. In the case of multiple treated observations
#' (you used \code{\link{scdataMulti}} to prepare the data), the user can specify \code{V} as a string equal to either "separate" or "pooled".
#' If \code{scdata()} was used to prepare the data, \code{V} is automatically set to "separate" as the two options are
#' equivalent. See the \strong{Details} section for more.
#' @param V.mat A conformable weighting matrix \eqn{\mathbf{V}} to be used in the minimization of the sum of squared residuals
#' \deqn{(\mathbf{A}-\mathbf{B}\mathbf{w}-\mathbf{C}\mathbf{r})'\mathbf{V}(\mathbf{A}-\mathbf{B}\mathbf{w}-\mathbf{C}\mathbf{r}).}
#' See the \strong{Details} section for more information on how to prepare this matrix.
#' @param solver a string containing the name of the solver used by \code{CVXR} when computing the weights. You can check which solvers are available
#' on your machine by running \code{CVXR::installed_solvers()}. More information on what different solvers do can be found
#' at the following link https://cvxr.rbind.io/cvxr_examples/cvxr_using-other-solvers/. "OSQP" is the default solver when 'lasso'
#' is the constraint type, whilst "ECOS" is the default in all other cases.
#' @param plot a logical specifying whether \code{\link{scplot}} should be called and a plot saved in the current working
#' directory. For more options see \code{\link{scplot}}.
#' @param plot.name a string containing the name of the plot (the format is by default .png). For more options see \code{\link{scplot}}.
#' @param plot.path a string containing the path at which the plot should be saved (default is output of \code{getwd()}.)
#' @param save.data a character specifying the name and the path of the saved dataframe containing the processed data used to produce the plot.
#'
#' @return
#' The function returns an object of class 'scest' containing two lists. The first list is labeled 'data' and
#' contains used data as returned by \code{\link{scdata}} and some other values.
#' \item{A}{a matrix containing pre-treatment features of the treated unit(s).}
#' \item{B}{a matrix containing pre-treatment features of the control units.}
#' \item{C}{a matrix containing covariates for adjustment.}
#' \item{P}{a matrix whose rows are the vectors used to predict the out-of-sample series for the synthetic unit(s).}
#' \item{P.diff}{for internal use only.}
#' \item{Y.pre}{a matrix containing the (raw) pre-treatment outcome of the treated unit(s).}
#' \item{Y.post}{a matrix containing the (raw) post-treatment outcome of the treated unit(s).}
#' \item{Y.pre.agg}{a matrix containing the aggregate pre-treatment outcome of the treated unit(s). This differs from
#' Y.pre only in the case 'effect' in \code{scdataMulti()} is set to either 'unit' or 'time'.}
#' \item{Y.post.agg}{a matrix containing the aggregate post-treatment outcome of the treated unit(s). This differs from
#' Y.post only in the case 'effect' in \code{scdataMulti()} is set to either 'unit' or 'time'.}
#' \item{Y.donors}{a matrix containing the pre-treatment outcome of the control units.}
#' \item{specs}{a list containing some specifics of the data:
#' \itemize{
#' \item{\code{J}, the number of control units}
#' \item{\code{K}, a numeric vector with the number of covariates used for adjustment for each feature}
#' \item{\code{M}, number of features}
#' \item{\code{KM}, the total number of covariates used for adjustment}
#' \item{\code{KMI}, the total number of covariates used for adjustment}
#' \item{\code{I}, number of treated units}
#' \item{\code{period.pre}, a numeric vector with the pre-treatment period}
#' \item{\code{period.post}, a numeric vector with the post-treatment period}
#' \item{\code{T0.features}, a numeric vector with the number of periods used in estimation for each feature}
#' \item{\code{T1.outcome}, the number of post-treatment periods}
#' \item{\code{constant}, for internal use only}
#' \item{\code{effect}, for internal use only}
#' \item{\code{anticipation}, number of periods of potential anticipation effects}
#' \item{\code{out.in.features}, for internal use only}
#' \item{\code{treated.units}, list containing the IDs of all treated units}
#' \item{\code{donors.list}, list containing the IDs of the donors of each treated unit}
#' \item{\code{class.type}, for internal use only}}}
#'
#' The second list is labeled 'est.results' and contains estimation results.
#' \item{w}{a matrix containing the estimated weights of the donors.}
#' \item{r}{a matrix containing the values of the covariates used for adjustment.}
#' \item{b}{a matrix containing \eqn{\mathbf{w}} and \eqn{\mathbf{r}}.}
#' \item{Y.pre.fit}{a matrix containing the estimated pre-treatment outcome of the SC unit(s).}
#' \item{Y.post.fit}{a matrix containing the estimated post-treatment outcome of the SC unit(s).}
#' \item{A.hat}{a matrix containing the predicted values of the features of the treated unit(s).}
#' \item{res}{a matrix containing the residuals \eqn{\mathbf{A}-\widehat{\mathbf{A}}}.}
#' \item{V}{a matrix containing the weighting matrix used in estimation.}
#' \item{w.constr}{a list containing the specifics of the constraint set used on the weights.}
#'
#' @details
#' Information is provided for the simple case in which \eqn{N_1=1} if not specified otherwise.
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
#'
#' \item{If \code{name == "L1-L2"}, then
#' \deqn{||\mathbf{w}||_1 = 1,\:\:\: ||\mathbf{w}||_2 \leq Q, \:\:\: w_j \geq 0,\:\: j =1,\ldots,J.}
#' where \eqn{Q} is a tuning parameter computed as in the "ridge" case.}
#' }}
#'
#' \item{\strong{Weighting Matrix.}
#' \itemize{
#' \item{if \code{V <- "separate"}, then \eqn{\mathbf{V} = \mathbf{I}} and the minimized objective function is
#' \deqn{\sum_{i=1}^{N_1} \sum_{l=1}^{M} \sum_{t=1}^{T_{0}}\left(a_{t, l}^{i}-\mathbf{b}_{t, l}^{{i \prime }} \mathbf{w}^{i}-\mathbf{c}_{t, l}^{{i \prime}} \mathbf{r}_{l}^{i}\right)^{2},}
#' which optimizes the separate fit for each treated unit.}
#' \item{if \code{V <- "pooled"}, then \eqn{\mathbf{V} = \frac{1}{I}\mathbf{1}\mathbf{1}'\otimes \mathbf{I}} and the minimized objective function is
#' \deqn{\sum_{l=1}^{M} \sum_{t=1}^{T_{0}}\left(\frac{1}{N_1^2} \sum_{i=1}^{N_1}\left(a_{t, l}^{i}-\mathbf{b}_{t, l}^{i \prime} \mathbf{w}^{i}-\mathbf{c}_{t, l}^{i\prime} \mathbf{r}_{l}^{i}\right)\right)^{2},}
#' which optimizes the pooled fit for the average of the treated units.}
#' \item{if the user wants to provide their own weighting matrix, then it must use the option \code{V.mat} to input a \eqn{v\times v} positive-definite matrix, where \eqn{v} is the 
#' number of rows of \eqn{\mathbf{B}} (or \eqn{\mathbf{C}}) after potential missing values have been removed. In case the user
#' wants to provide their own \code{V}, we suggest to check the appropriate dimension \eqn{v} by inspecting the output
#' of either \code{scdata} or \code{scdataMulti} and check the dimensions of \eqn{\mathbf{B}} (and \eqn{\mathbf{C}}). Note that
#' the weighting matrix could cause problems to the optimizer if not properly scaled. For example, if \eqn{\mathbf{V}} is diagonal
#' we suggest to divide each of its entries by \eqn{\|\mathrm{diag}(\mathbf{V})\|_1}.}
#' }}
#' }
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
#' @seealso \code{\link{scdataMulti}}, \code{\link{scdata}}, \code{\link{scpi}}, \code{\link{scplot}}, \code{\link{scplotMulti}}
#'
#' @examples
#'
#' data <- scpi_germany
#'
#' df <- scdata(df = data, id.var = "country", time.var = "year",
#'              outcome.var = "gdp", period.pre = (1960:1990),
#'              period.post = (1991:2003), unit.tr = "West Germany",
#'              unit.co = setdiff(unique(data$country), "West Germany"),
#'              constant = TRUE, cointegrated.data = TRUE)
#'
#' result <- scest(df, w.constr = list(name = "simplex", Q = 1))
#' result <- scest(df, w.constr = list(lb = 0, dir = "==", p = "L1", Q = 1))
#'
#' @export

scest <- function(data,
                  w.constr  = NULL,
                  V         = "separate",
                  V.mat     = NULL,
                  solver    = "ECOS",
                  plot      = FALSE,
                  plot.name = NULL,
                  plot.path = NULL,
                  save.data = NULL) {

  ##########################
  if ((methods::is(data, "scdata") || methods::is(data, "scdataMulti")) == FALSE) {
    stop("data should be the object returned by running scdata or scdata_multi!")
  }

  if (is.null(w.constr) == FALSE) { # The user has specified W

    if (is.list(w.constr) == FALSE) {
      stop("w.constr should be a list!")
    }

    if (!"name" %in% names(w.constr)) {
      w.constr[["name"]] <- "NULL"
    } else {
      if (!w.constr[["name"]] %in% c("simplex","lasso","ridge","ols","L1-L2")) {
        stop("If 'name' is specified in w.constr, then it should be chosen among
             'simplex', 'lasso', 'ridge', 'ols', 'L1-L2'.")
      }
    }
  }

  if (!(solver %in% CVXR::installed_solvers())) {
    stop(paste0("The specified solver - ", solver," - is not available on your machine! Run
                CVXR::installed_solvers() to see the list of available options on your machine."))
  }

  # Data matrices
  A <- data$A
  B <- data$B
  C <- data$C
  P <- data$P
  Z <- cbind(B, C)
  Y.donors <- data$Y.donors
  outcome.var <- data$specs$outcome.var
  
  if (class(data)[1] == 'scdata') {
    class.type <- 'scpi_data'
  } else if (class(data)[1] == 'scdataMulti') {
    class.type <- 'scpi_data_multi'
  }
  
  V.type <- V
  
  # Data specs
  if (class.type == 'scpi_data') {
    J               <- data$specs$J
    KM              <- data$specs$KM
    I               <- 1
    KMI             <- KM
    Jtot            <- J
    M               <- data$specs$M
    
  } else if (class.type == 'scpi_data_multi') {
    J               <- data$specs$J # total number of donors
    Jtot            <- sum(unlist(J))
    KM              <- data$specs$KM # total number of covs used for adjustment per treated unit
    KMI             <- data$specs$KMI # total number of covariates used for adjustment
    I               <- data$specs$I  # number of treated units
    T0.M            <- lapply(data$specs$T0.features, sum) # observations per treated unit
    M               <- data$specs$M
  }

  T0.features     <- data$specs$T0.features
  out.in.features <- data$specs$out.in.features

  ##########################
  ## Set up the estimation problem

  # Create weighting matrix
  if (is.character(V) == FALSE) {
    stop("The object V should be a string! If you want to specify manually the weighting matrix use the option V.mat!")
  }
  
  if (is.null(V.mat) == FALSE) {
    if (is.matrix(V.mat) == FALSE) {
      stop("The object V.mat should be a matrix!")
    }

    if (nrow(V.mat) != nrow(B) || ncol(V.mat) != nrow(B)) {
      stop(paste0("V.mat should be a ", nrow(B), "x", nrow(B)," matrix, but currently it is
                a ", nrow(V.mat), "x", ncol(V.mat), " matrix!"))
    }
    rownames(V.mat) <- rownames(Z)
    colnames(V.mat) <- rownames(V.mat)
    
  } else {
    V.mat <- V.prep(type = V.type, B, T0.features, I)
  }

  V <- V.mat

  # Create lists of matrices
  # Estimate SC
  if (class.type == 'scpi_data') { # single treated unit
    w.constr <- w.constr.OBJ(w.constr, A, Z, V, J, KM, M)
    if (w.constr[["name"]] == "lasso") solver <- "OSQP"
    b <- b.est(A = A, Z = Z, J = J, KM = KM, w.constr = w.constr, V = V, CVXR.solver = solver)

  } else if (class.type == 'scpi_data_multi') { # multiple treated units
    
    A.list <- mat2list(A)
    B.list <- mat2list(B)
    C.list <- mat2list(C)
    V.list <- mat2list(V)

    w.constr.list <- list()

    w.store <- c()
    r.store <- c()
    j.lb <- 1
    j.ub <- J[[1]]

    for (i in seq_len(I)) {
      if (i > 1){
        j.lb <- j.ub + 1
        j.ub <- j.lb + J[[i]] - 1
      }

      A.i <- A.list[[i]]
      Z.i <- cbind(B.list[[i]], C.list[[i]])
      V.i <- V.list[[i]]

      w.constr.list[[data$specs$treated.units[i]]] <- w.constr.OBJ(w.constr, A.i, Z.i, V.i, J[[i]], KM[[i]], M[[i]])
      if (w.constr.list[[i]]["name"] == "lasso") solver <- "OSQP"

      if (V.type == "separate") {
        res <- b.est(A = A.i, Z = Z.i, J = J[[i]], KM = KM[[i]], w.constr = w.constr.list[[i]], V = V.i, CVXR.solver = solver)
        w.store <- c(w.store, res[1:J[[i]]])
        if (KM[[i]] > 0) r.store <- c(r.store, res[(J[[i]]+1):length(res)])
      }
    }

    if (V.type != "separate") {
      if (w.constr.list[[i]]["name"] == "lasso") solver <- "OSQP"
      b <- b.est.multi(A = A, Z = Z, J = J, KMI = KMI, I = I, 
                       w.constr = w.constr.list, V = V, CVXR.solver = solver)
      b <- b[,1,drop=TRUE]
      
    } else if (V.type == "separate") {
      b <- c(w.store, r.store)
    }
    w.constr <- w.constr.list
  }

  ##########################
  ## Create useful objects
  # Store point estimates
  if (KMI == 0) {
    w <- b                       # Estimated weights
    r <- NULL                    # Loading factors of additional covariates
  } else {
    w <- b[1:Jtot]
    r <- b[(Jtot + 1):length(b)]
  }

  # Fitted values and residuals
  A.hat    <- Z %*% b                            # Pre-treatment fit of synthetic unit
  res      <- A  -  A.hat                        # Estimation residual u

  # Pre-treatment fit of outcome of interest
  if (class.type == 'scpi_data'){
    if (out.in.features == TRUE) {
      fit.pre  <- A.hat[1:T0.features[outcome.var], , drop = FALSE]
      names    <- strsplit(rownames(fit.pre), "\\.")
      rownames(fit.pre) <- unlist(lapply(names, function(x) paste(x[1],x[3],sep=".")))
    } else {
      fit.pre  <- Y.donors %*% w
    }
  } else if(class.type == 'scpi_data_multi') {
    
    i.lb <- 1
    fit.pre <- c()
    Yd.list <- mat2list(data$Y.donors)
    w.list <- mat2list(as.matrix(w), cols=FALSE)
    
    for (i in seq_len(I)) {
      i.ub <- i.lb + T0.features[[i]][outcome.var] - 1
      if (out.in.features[[i]] == TRUE) {
        fit.pre.i  <- A.hat[i.lb:i.ub, , drop = FALSE]
        names    <- strsplit(rownames(fit.pre.i), "\\.")
        rownames(fit.pre.i) <- unlist(lapply(names, function(x) paste(x[1],x[3],sep=".")))   
      } else {
        fit.pre.i <- Yd.list[[i]] %*% w.list[[i]]
      }
      i.lb <- i.lb + sum(unlist(T0.features[[i]]), na.rm = TRUE)
      fit.pre <- rbind(fit.pre, fit.pre.i)
    }
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

  if (class.type == 'scpi_data') {
    df   <- list(A = data$A,
                 B = data$B,
                 C = data$C,
                 P = data$P,
                 Z = Z,
                 specs  = data$specs,
                 Y.pre  = data$Y.pre,
                 Y.post = data$Y.post,
                 Y.pre.agg = data$Y.pre,
                 Y.post.agg = data$Y.post.agg)

  } else if (class.type == 'scpi_data_multi') {
    
    # shortcut to avoid "no visible binding for global variable 'X' when checking the package
    Treatment <- NULL
    
    # Y.pre and Y.post might require some extra work
    # if the predictand of interest is aggregate (either over units or over time)
    # The next function process the data in the same way it's done in scplotMulti
    treated.units   <- data$specs$treated.units
    sparse.matrices <- data$specs$sparse.matrices
    anticipation    <- data$specs$anticipation
    period.post     <- data$specs$period.post
    units.est       <- data$specs$units.est
    effect          <- data$specs$effect
    Y.df            <- data$Y.df
    Y.pre.fit       <- fit.pre
    Y.post.fit      <- fit.post
    
    # create to plot object
    Yprocessed <- outcomeGet(Y.pre.fit=Y.pre.fit, Y.post.fit=Y.post.fit, Y.df=Y.df,
                             units.est=units.est, treated.units=treated.units, plot.type=effect,
                             anticipation=anticipation, period.post=period.post,
                             sparse.matrices=sparse.matrices)

    Ydf.pre <- subset(Yprocessed$toplot, Treatment == 0)
    Ydf.post <- subset(Yprocessed$toplot, Treatment == 1)
    Y.pre.agg  <- as.matrix(Ydf.pre[["Actual"]])
    Y.post.agg <- as.matrix(Ydf.post[["Actual"]])
    names.pre <- paste(Ydf.pre$ID, Ydf.pre$Time, sep=".")
    names.post <- paste(Ydf.post$ID, Ydf.post$Time, sep=".")
    rownames(Y.pre.agg) <- names.pre
    rownames(Y.post.agg) <- names.post
    colnames(Y.pre.agg) <- data$specs$outcome.var
    colnames(Y.post.agg) <- data$specs$outcome.var
    
    df   <- list(A = data$A,
                 B = data$B,
                 C = data$C,
                 P = data$P,
                 P.diff = data$P.diff,
                 Y.df = data$Y.df,
                 Y.pre = data$Y.pre,
                 Y.post = data$Y.post,
                 Y.pre.agg = Y.pre.agg,
                 Y.post.agg = Y.post.agg,
                 Z = Z,
                 specs  = data$specs)
  }

  ##########################
  ## Return to the user
  to.return <- list(data = df, est.results = est.results)
  class(to.return) <- 'scest'
  if (class.type == 'scpi_data') {
    to.return$data$specs$class.type <- 'scpi_scest'
  } else if (class.type == 'scpi_data_multi') {
    to.return$data$specs$class.type <- 'scpi_scest_multi'
  }

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

    if (class.type == 'scpi_data'){

      scplot(result = to.return, fig.path = fig.path,
             fig.name = fig.name, fig.format = "png", save.data = save.data)

      } else if (class.type == "scpi_data_multi"){

      scplotMulti(result = to.return)

    }

  }

  return(to.return)
}

###############################################################################
