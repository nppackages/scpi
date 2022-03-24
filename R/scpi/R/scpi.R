###############################################################################

#' @title Prediction Intervals for Synthetic Control Methods
#'
#' @description The command implements estimation and inference procedures for Synthetic Control (SC) methods using least square, lasso, ridge, or simplex-type constraints. Uncertainty is quantified using prediction
#' intervals according to \href{https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, M. D., Feng, Y., & Titiunik, R. (2021)}. \code{\link{scpi}} returns the estimated 
#' post-treatment series for the synthetic unit through the command \code{\link{scest}} and quantifies in-sample and out-of-sample uncertainty to provide confidence intervals
#' for each point estimate.
#' 
#' Companion \href{https://www.stata.com/}{Stata} and \href{https://www.python.org/}{Python} packages are described in \href{https://arxiv.org/abs/2202.05984}{Cattaneo, Feng, Palomba, and Titiunik (2022)}.
#'
#' Companion commands are: \link{scdata} for data preparation, \link{scest} for point estimation, and \link{scplot} for plots.
#' 
#' Related Stata, R, and Python packages useful for inference in SC designs are described in the following website:
#' 
#' \href{ https://nppackages.github.io/scpi/}{ https://nppackages.github.io/scpi/}
#' 
#' For an introduction to synthetic control methods, see \href{https://www.aeaweb.org/articles?id=10.1257/jel.20191450}{Abadie (2021)} and references therein.
#'
#' @param data a class `scpi_data' object, obtained by calling \code{\link{scdata}}.
#'
#' @param w.constr a list specifying the constraint set the estimated weights of the donors must belong to.
#' \code{w.constr} can contain up to five elements:
#' - `\code{p}', a scalar indicating the norm to be used (\code{p} should be one of "no norm", "L1", and "L2")
#' - `\code{dir}', a string indicating whether the constraint on the norm is an equality ("==") or inequality ("<=")
#' - `\code{Q}', a scalar defining the value of the constraint on the norm
#' - `\code{lb}', a scalar defining the lower bound on the weights. It can be either 0 or \code{-Inf}.
#' - `\code{name}', a character selecting one of the default proposals
#' See the \strong{Details} section for more.
#'
#' @param V a weighting matrix to be used when minimizing
#' \deqn{(\mathbf{A}-\mathbf{B}\mathbf{w}-\mathbf{C}\mathbf{r})'\mathbf{V}(\mathbf{A}-\mathbf{B}\mathbf{w}-\mathbf{C}\mathbf{r})}
#' @param P a \eqn{T_1\times (J+K_1)} matrix containing the design matrix to be used to obtain the predicted.
#' post-intervention outcome of the synthetic control unit. \eqn{T_1} is the number of post-treatment periods,
#' \eqn{J} is the size of the donor pool, and \eqn{K_1} is the number of covariates used for adjustment in the outcome equation.
#' @param rho a string specifying the regularizing parameter that imposes sparsity on the estimated vector of weights. If
#' \code{rho = 'type-1'} (the default), then the tuning parameter is computed based on optimization inequalities. Users can provide a scalar 
#' with their own value for \code{rho}. Other options are described in the \strong{Details} section.
#' @param rho.max a scalar indicating the maximum value attainable by the tuning parameter \code{rho}.
#' @param u.missp a logical indicating if misspecification should be taken into account when dealing with \eqn{\mathbf{u}}.
#' @param u.order a scalar that sets the order of the polynomial in \eqn{\mathbf{B}} when predicting moments of \eqn{\mathbf{u}}. 
#' The default is \code{u.order = 1}, however if there is risk of over-fitting, the command automatically sets it 
#' to \code{u.order = 0}. See the \strong{Details} section for more information.
#' @param u.lags a scalar that sets the number of lags of \eqn{\mathbf{B}} when predicting moments of \eqn{\mathbf{u}}.
#' The default is \code{u.lags = 0}, however if there is risk of over-fitting, the command automatically sets it 
#' to \code{u.lags = 0}. See the \strong{Details} section for more information.
#' @param u.design a matrix with the same number of rows of \eqn{\mathbf{A}} and \eqn{\mathbf{B}} and whose columns specify the design matrix
#' to be used when modeling the estimated pseudo-true residuals \eqn{\mathbf{u}}.
#' @param u.sigma a string specifying the type of variance-covariance estimator to be used when estimating
#'    the conditional variance of \eqn{\mathbf{u}}.
#' @param u.alpha a scalar specifying the confidence level for in-sample uncertainty, i.e. 1 - \code{u.alpha} is the confidence level.
#' @param e.method a string selecting the method to be used in quantifying out-of-sample uncertainty among:
#'  "gaussian" which uses conditional subgaussian bounds; "ls" which specifies a location-scale model for \eqn{\mathbf{u}}; "qreg" which employs a
#'  quantile regressions to get the conditional bounds; "all" uses each one of the previous methods.
#' @param e.order a scalar that sets the order of the polynomial in \eqn{\mathbf{B}} when predicting moments of \eqn{\mathbf{e}}.
#' The default is \code{e.order = 1}, however if there is risk of over-fitting, the command automatically sets it 
#' to \code{e.order = 0}. See the \strong{Details} section for more information.
#' @param e.lags a scalar that sets the number of lags of \eqn{\mathbf{B}} when predicting moments of \eqn{\mathbf{e}}.
#' The default is \code{e.order = 1}, however if there is risk of over-fitting, the command automatically sets it 
#' to \code{e.order = 0}. See the \strong{Details} section for more information.
#' @param e.design a matrix with the same number of rows of \eqn{\mathbf{A}} and \eqn{\mathbf{B}} and whose columns specify the design matrix
#' to be used when modeling the estimated out-of-sample residuals \eqn{\mathbf{e}}.
#' @param e.alpha a scalar specifying the confidence level for out-of-sample uncertainty, i.e. 1 - \code{e.alpha} is the confidence level.
#' @param sims a scalar providing the number of simulations to be used in quantifying in-sample uncertainty.
#' @param plot a logical specifying whether \code{\link{scplot}} should be called and a plot saved in the current working directory. For more options see \code{\link{scplot}}.
#' @param plot.name a string containing the name of the plot (the format is by default .png). For more options see \code{\link{scplot}}.
#' @param cores number of cores to be used by the command. The default is the number of cores available minus one.
#' @param w.bounds a \eqn{T_1\times 2} matrix with the user-provided bounds on \eqn{\beta}. If \code{w.bounds} is provided, then
#' the quantification of in-sample uncertainty is skipped. It is possible to provide only the lower bound or the upper bound
#' by filling the other column with \code{NA}s.
#' @param e.bounds a \eqn{T_1\times 2} matrix with the user-provided bounds on \eqn{(\widehat{\mathbf{w}},\widehat{\mathbf{r}})^{\prime}}. If \code{e.bounds} is provided, then
#' the quantification of out-of-sample uncertainty is skipped. It is possible to provide only the lower bound or the upper bound
#' by filling the other column with \code{NA}s.
#' @param opt.list.est a list specifying the stopping criteria and the algorithm for the underlying optimizer \code{\link{nloptr}} or \code{\link{CVXR}} for point estimation.  
#' See the \strong{Details} section for more. 
#' @param opt.list.inf similar to the previous one but for the optimizer used for inference purposes. See the \strong{Details} section for more. 
#'  
#' @param save.data a character specifying the name and the path of the saved dataframe containing the processed data used to produce the plot. 
#' 
#' @param verbose if \code{TRUE} prints additional information in the console.
#' 
#' @return
#' The function returns an object of class 'scpi_scpi' containing three lists. The first list is labeled 'data' and contains used data as returned by \code{\link{scdata}} and some other values.
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
#' The second list is labeled 'est.results' containing all the results from \code{\link{scest}}.
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
#' The third list is labeled 'inference.results' and contains all the inference-related results.
#' \item{CI.in.sample}{a matrix containing the prediction intervals taking only in-sample uncertainty in to account.}
#' \item{CI.all.gaussian}{a matrix containing the prediction intervals estimating out-of-sample uncertainty with sub-Gaussian bounds.}
#' \item{CI.all.ls}{a matrix containing the prediction intervals estimating out-of-sample uncertainty with a location-scale model.}
#' \item{CI.all.qreg}{a matrix containing the prediction intervals estimating out-of-sample uncertainty with quantile regressions.}
#' \item{Sigma}{a matrix containing the estimated (conditional) variance-covariance \eqn{\boldsymbol{\Sigma}}.}
#' \item{u.mean}{a matrix containing the estimated (conditional) mean of the pseudo-residuals \eqn{\mathbf{u}}.}
#' \item{u.var}{a matrix containing the estimated (conditional) variance-covariance of the pseudo-residuals \eqn{\mathbf{u}}.}
#' \item{e.mean}{a matrix containing the estimated (conditional) mean of the out-of-sample error \eqn{e}.}
#' \item{e.var}{a matrix containing the estimated (conditional) variance of the out-of-sample error \eqn{e}.}
#' \item{u.missp}{a logical indicating whether the model has been treated as misspecified or not.}
#' \item{u.lags}{an integer containing the number of lags in B used in predicting moments of the pseudo-residuals \eqn{\mathbf{u}}.}
#' \item{u.order}{an integer containing the order of the polynomial in B used in predicting moments of the pseudo-residuals \eqn{\mathbf{u}}.}
#' \item{u.sigma}{a string indicating the estimator used for \code{Sigma}.}
#' \item{u.user}{a logical indicating whether the design matrix to predict moments of \eqn{\mathbf{u}} was user-provided.}
#' \item{u.T}{a scalar indicating the number of observations used to predict moments of \eqn{\mathbf{u}}.}
#' \item{u.params}{a scalar indicating the number of parameters used to predict moments of \eqn{\mathbf{u}}.}
#' \item{u.D}{the design matrix used to predict moments of \eqn{\mathbf{u}},}
#' \item{u.alpha}{a scalar determining the confidence level used for in-sample uncertainty, i.e. 1-\code{u.alpha} is the confidence level.}
#' \item{e.method}{a string indicating the specification used to predict moments of the out-of-sample error \eqn{e}.}
#' \item{e.lags}{an integer containing the number of lags in B used in predicting moments of the out-of-sample error \eqn{e}.}
#' \item{e.order}{an integer containing the order of the polynomial in B used in predicting moments of the out-of-sample error \eqn{e}.}
#' \item{e.user}{a logical indicating whether the design matrix to predict moments of \eqn{e} was user-provided.}
#' \item{e.T}{a scalar indicating the number of observations used to predict moments of \eqn{\mathbf{u}}.}
#' \item{e.params}{a scalar indicating the number of parameters used to predict moments of \eqn{\mathbf{u}}.}
#' \item{e.alpha}{a scalar determining the confidence level used for out-of-sample uncertainty, i.e. 1-\code{e.alpha} is the confidence level.}
#' \item{e.D}{the design matrix used to predict moments of \eqn{\mathbf{u}},}
#' \item{rho}{an integer specifying the estimated regularizing parameter that imposes sparsity on the estimated vector of weights.}
#' \item{sims}{an integer indicating the number of simulations used in quantifying in-sample uncertainty.}
#' \item{failed.sims}{a matrix containing the number of failed simulations per post-treatment period to estimate lower and upper bounds.}
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
#' \item{\strong{Regularization.} \code{rho} is estimated through the formula
#' \deqn{\varrho = \mathcal{C}\frac{\log (T_0)^c}{T_0^{1/2}}}
#' where \eqn{\mathcal{C} = \widehat{\sigma}_u / \min_j \widehat{\sigma}_{b_j}} if \code{rho = 'type-1'}, 
#' \eqn{\mathcal{C} = \max_{j}\widehat{\sigma}_{b_j}\widehat{\sigma}_{u} / \min_j \widehat{\sigma}_{b_j}^2} if \code{rho = 'type-2'}, and
#' \eqn{\mathcal{C} = \max_{j}\widehat{\sigma}_{b_ju} / \min_j \widehat{\sigma}_{b_j}^2} if \code{rho = 'type-3'},
#' 
#' \code{rho} defines a new sparse weight vector as
#'  \deqn{\widehat{w}^\star_j = \mathbf{1}(\widehat{w}_j\geq \varrho)}
#' }
#'
#' \item{\strong{In-sample uncertainty.} To quantify in-sample uncertainty it is necessary to model the pseudo-residuals \eqn{\mathbf{u}}.
#' First of all, estimation of the first moment of \eqn{\mathbf{u}} can be controlled through
#' the option \code{u.missp}. When \code{u.missp = FALSE}, then \eqn{\mathbf{E}[u\: |\: \mathbf{D}_u]=0}. If instead \code{u.missp = TRUE},
#' then \eqn{\mathbf{E}[u\: |\: \mathbf{D}_u]} is estimated using a linear regression of
#' \eqn{\widehat{\mathbf{u}}} on \eqn{\mathbf{D}_u}. The default set of variables in \eqn{\mathbf{D}_u} is composed of \eqn{\mathbf{B}}, 
#' \eqn{\mathbf{C}} and, if required, it is augmented with lags (\code{u.lags}) and polynomials (\code{u.order}) of \eqn{\mathbf{B}}. 
#' The option \code{u.design} allows the user to provide an ad-hoc set of variables to form \eqn{\mathbf{D}_u}. 
#' Regarding the second moment of \eqn{\mathbf{u}}, different estimators can be chosen:
#' HC0, HC1, HC2, HC3, and HC4 using the option \code{u.sigma}.}
#'
#' \item{\strong{Out-of-sample uncertainty.} To quantify out-of-sample uncertainty it is necessary to model the out-of-sample residuals
#' \eqn{\mathbf{e}} and estimate relevant moments. By default, the design matrix used during estimation \eqn{\mathbf{D}_e} is composed of the blocks in 
#' \eqn{\mathbf{B}} and \eqn{\mathbf{C}} corresponding to the outcome variable. Moreover, if required by the user, \eqn{\mathbf{D}_e}
#' is augmented with lags (\code{e.lags}) and polynomials (\code{e.order}) of \eqn{\mathbf{B}}. The option \code{e.design} allows the user to provide an
#' ad-hoc set of variables to form \eqn{\mathbf{D}_e}. Finally, the option \code{e.method} allows the user to select one of three
#' estimation methods: "gaussian" relies on conditional sub-Gaussian bounds; "ls" estimates conditional bounds using a location-scale
#' model; "qreg" uses conditional quantile regression of the residuals \eqn{\mathbf{e}} on \eqn{\mathbf{D}_e}.}
#' 
#' \item{\strong{Residual Estimation Over-fitting.} To estimate conditional moments of \eqn{\mathbf{u}} and \eqn{e_t} 
#' we rely on two design matrices, \eqn{\mathbf{D}_u} and \eqn{\mathbf{D}_e} (see above). Let \eqn{d_u} and \eqn{d_e} be the number of 
#' columns in \eqn{\mathbf{D}_u} and \eqn{\mathbf{D}_e}, respectively. Assuming no missing values and balanced features, the
#' number of observation used to estimate moments of \eqn{\mathbf{u}} is \eqn{T_0\cdot M}, whilst for moments of \eqn{e_t} is \eqn{T_0}.
#' Our rule of thumb to avoid over-fitting is to check if \eqn{T_0\cdot M \geq d_u + 20} or \eqn{T_0 \geq d_e + 20}. If the 
#' former condition is not satisfied we automatically set \code{u.order = u.lags = 0}, if instead the latter is not met 
#' we automatically set \code{e.order = e.lags = 0}.}
#' 
#' \item{\strong{Algorithm Options.} The default is a sequential quadratic programming (SQP) algorithm for nonlinearly constrained gradient-based optimization 
#' (\code{algorithm = 'NLOPTR_LD_SLSQP'}) for all cases not involving the L1 norm. 
#' For a complete list of algorithms see \href{https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/}{the official \code{nloptr} documentation}.
#' The other default values are \code{maxeval = 5000}, \code{ftol_res = 1.0e-8}, \code{ftol_abs = 1.0e-8}, \code{xtol_res = 1.0e-8}, \code{xtol_abs = 1.0e-8},
#' \code{tol_constraints_eq = 1.0e-8}, and, finally, \code{tol_constraints_ineq = 1.0e-8}. More information on the stopping criteria can be obtained running
#' \code{nloptr.print.options()} or \code{nloptr.get.default.options()}. If the optimization involves the L1 norm then \code{CVXR} is used for optimization.  
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
#' \item{\href{https://www.aeaweb.org/articles?id=10.1257/jel.20191450}{Abadie, A. (2021)}. Using synthetic controls: Feasibility, data requirements, and methodological aspects.
#' \emph{Journal of Economic Literature}, 59(2), 391-425.}
#' \item{\href{https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, M. D., Feng, Y., & Titiunik, R. 
#' (2021)}. Prediction intervals for synthetic control methods. \emph{Journal of the American Statistical Association}, 116(536), 1865-1880.}
#' \item{\href{https://arxiv.org/abs/2202.05984}{Cattaneo, M. D., Feng, Y., Palomba F., and Titiunik, R. (2022).}
#' scpi: Uncertainty Quantification for Synthetic Control Estimators, \emph{arXiv}:2202.05984.}
#' }
#'
#' @seealso \code{\link{scdata}}, \code{\link{scest}}, \code{\link{scplot}}
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
#' result <- scpi(df, w.constr = list(name = "simplex", Q = 1), cores = 1, sims = 10)
#' result <- scpi(df, w.constr = list(lb = 0, dir = "==", p = "L1", Q = 1),
#'                cores = 1, sims = 10)
#'                            
#' @export

scpi  <- function(data,
                  w.constr     = NULL,
                  V            = NULL,
                  P            = NULL,
                  u.missp      = TRUE,
                  u.sigma      = "HC1",
                  u.order      = 1,
                  u.lags       = 0,
                  u.design     = NULL,
                  u.alpha      = 0.05,
                  e.method     = "all",
                  e.order      = 1,
                  e.lags       = 0,
                  e.design     = NULL,
                  e.alpha      = 0.05,
                  sims         = 200,
                  rho          = NULL,
                  rho.max      = NULL,
                  cores        = NULL,
                  plot         = FALSE,
                  plot.name    = NULL,
                  w.bounds     = NULL,
                  e.bounds     = NULL,
                  opt.list.est = NULL,
                  opt.list.inf = NULL,
                  save.data    = NULL,
                  verbose      = TRUE) {


  if (methods::is(data, "scpi_data") == FALSE) {
    stop("data should be the object returned by running scdata!!")
  }

  #############################################################################
  #############################################################################
  ## Estimation of synthetic weights
  if (verbose) {
    cat("---------------------------------------------------------------\n")
    cat("Estimating Weights...\n")
  }
  sc.pred <- scest(data = data, w.constr = w.constr, V = V, opt.list = opt.list.est)


  #############################################################################
  #############################################################################
  ## Retrieve processed data from scest

  A          <- sc.pred$data$A                           # Features of treated unit
  B          <- sc.pred$data$B                           # Features of control units
  C          <- sc.pred$data$C                           # Covariates for adjustment
  Z          <- sc.pred$data$Z                           # B and C column-bind
  Y.donors   <- data$Y.donors                            # Outcome variable of control units
  K          <- sc.pred$data$specs$K                     # Number of covs for adjustment per feature
  KM         <- sc.pred$data$specs$KM                    # Dimension of r (total number of covs for adj)
  J          <- sc.pred$data$specs$J                     # Number of donors
  M          <- sc.pred$data$specs$M                     # Number of features
  T0         <- sc.pred$data$specs$T0.features           # Time periods used per feature
  T1         <- sc.pred$data$specs$T1.outcome            # Number of out-of-sample periods
  features   <- sc.pred$data$specs$features              # Name of features
  constant   <- sc.pred$data$specs$glob.cons             # Logical indicating whether a constant is included
  out.feat   <- sc.pred$data$specs$out.in.features       # Logical indicating whether the outcome variable is among features
  coig.data  <- sc.pred$data$specs$cointegrated.data     # Logical indicating whether B is cointegrated
  w.constr   <- sc.pred$est.results$w.constr             # Constraints on w
  V          <- sc.pred$est.results$V                    # Weighting matrix
  w          <- sc.pred$est.results$w                    # Estimated vector of weights
  r          <- sc.pred$est.results$r                    # Estimated coefficients of covariates
  b          <- sc.pred$est.results$b                    # w and r column-bind
  Y.post.fit <- sc.pred$est.results$Y.post.fit           # Estimated post-treatment outcome for SC unit
  res        <- sc.pred$est.results$res                  # Residuals from estimation
  T0.tot     <- sum(T0)                                  # Total number of observations used in estimation

  
  # Prepare algorithm options for inference (the one for estimation is handled in scest)
  if (is.null(opt.list.inf) == FALSE) {
    if (is.list(opt.list.inf) == FALSE) {
      stop("The object opt.list.inf should be a list!")
    }
  }

  if (is.null(opt.list.est) == FALSE) {
    if (is.list(opt.list.est) == FALSE) {
      stop("The object opt.list.est should be a list!")
    }
  }
  
  # Check on P
  if (is.null(P) == TRUE) {
    P <- sc.pred$data$P                         # Matrix for out-of-sample prediction

  } else { # User-provided prediction matrix P (should be T1 by (J+KM))
    if (is.matrix(P) == FALSE) {
      stop("The object P should be a matrix!")
    }
    if (nrow(P) != T1) {
      stop(paste("The matrix P currently has",nrow(P),"rows when instead",T1,"where expected
                 (i.e. the number of post-intervention periods)!"))
    }
    if (ncol(P) != (J+K[1])) {
      stop(paste("The matrix P currently has",ncol(P),"columns when instead",(J+K[1]),"where expected
                 (i.e. the size of the donor pool plus the number of covariates used in adjustment in the outcome equation)!"))
    }
    # Add zeros to avoid loading the coefficient of covs used for adj in other eqs
    if ((KM-K[1]) > 0) {
      zeros <- matrix(0, nrow(P), (KM-K[1]))
      P <- cbind(P, zeros)
    }
  }


  # Check on u_sigma and user-provided bounds
  if (is.character(u.sigma) == FALSE) {
    stop("The object u.sigma should be of type character!!")

  } else {
    if (!(u.sigma %in% c('HC0','HC1','HC2','HC3','HC4')))  {
      stop("Supported variance estimators are 'HC0','HC1','HC2','HC3','HC4'.")
    }
  } 
  
  if (is.null(w.bounds) == FALSE) {
    if (is.matrix(w.bounds) == FALSE) {
      stop("The object w.bounds should be a matrix!")
    }
    
    if (ncol(w.bounds) != 2) {
      stop("w.bounds should be a matrix with two columns: the first column for the 
           lower bound, the second for the upper. In case you don't want to specify
           the lower or the upper bound just fill the specific column with NAs.")
    }
    
    if (nrow(w.bounds) != length(Y.post.fit)) {
      stop(paste("w.bounds should be a matrix with ", length(Y.post.fit), 
                 " rows (i.e. the number of post-intervention periods)."))
    }
  }
  
  if (sims < 10) {
    stop("The number of simulations needs to be larger or equal than 10!")
  }
  
  if (is.null(e.bounds) == FALSE) {
    if (is.matrix(e.bounds) == FALSE) {
      stop("The object e.bounds should be a matrix!")
    }
    
    if (ncol(e.bounds) != 2) {
      stop("e.bounds should be a matrix with two columns: the first column for the 
           lower bound, the second for the upper. In case you don't want to specify
           the lower or the upper bound just fill the specific column with NAs.")
    }
    
    if (nrow(e.bounds) != length(Y.post.fit)) {
      stop(paste("e.bounds should be a matrix with ", length(Y.post.fit), 
                 " rows (i.e. the number of post-intervention periods)."))
    }
  }  
  
  # Check rho
  if (is.null(rho) == FALSE) {
    if (is.character(rho) == TRUE) {
      if(!(rho %in% c('type-1','type-2','type-3'))) {
        stop("When not a scalar, 'rho' must be 'type-1', 'type-2', or 'type-3'.")
      }
    }
  } else {
    rho <- "type-1"
  }
  
  # Check on number of cores
  if (is.null(cores) == FALSE) {
    n.cores <- detectCores(logical = TRUE)
    if (cores > n.cores) {
      stop(paste("You selected",cores,"cores, but only",n.cores," cores are available on your machine!"))
    }
  } else {
    cores <- detectCores(logical = TRUE) - 1
    warning(paste("scpi is using",cores,"cores for estimation! You can adjust the number of cores with the 'cores' option."), immediate. = TRUE, call. = FALSE)
  }
  
  if (verbose) {
    cat("Quantifying Uncertainty\n")
    executionTime(T0, J, T1, sims, cores, w.constr[['name']])
  }
  
  #############################################################################
  #############################################################################
  ### Estimate In-Sample Uncertainty
  #############################################################################
  #############################################################################
  
  ## Regularize W and local geometry
  loc.geom <- local.geom(w.constr, rho, rho.max, res, B, C, coig.data, T0.tot, J, w, verbose)
  
  w.star       <- loc.geom$w.star
  index.w      <- loc.geom$index.w
  w.constr.inf <- loc.geom$w.constr
  rho          <- loc.geom$rho
  Q.star       <- loc.geom$Q.star

  # Create an index that selects all non-zero weights and additional covariates
  index <- c(index.w, rep(TRUE, KM))
  beta  <- c(w.star, r)
  
  #############################################################################
  ##########################################################################
  ## Prepare design matrix for in-sample uncertainty
  obj <- u.des.prep(B, C, u.order, u.lags, coig.data, T0.tot, M, constant,
                        index, index.w, features, u.design, res)
  u.des.0 <- obj$u.des.0
  f.id <- as.factor(obj$f.id)
  
  #############################################################################
  ##########################################################################
  ## Prepare design matrices for out-of-sample uncertainty
  e.des <- e.des.prep(B, C, P, e.order, e.lags, res, sc.pred, Y.donors, out.feat, features,
                      J, M, index, index.w, coig.data, T0, T0.tot, T1, constant, e.design)
  
  e.res   <- e.des$e.res
  e.des.0 <- e.des$e.des.0
  e.des.1 <- e.des$e.des.1
  
  #############################################################################
  ###########################################################################
  # Remove NA - In Sample Uncertainty
  X  <- cbind(A, res, u.des.0, Z, f.id)
  XX <- na.omit(X)
  j1 <- 1
  j2 <- 2
  j3 <- j2 + 1
  j4 <- j2 + ncol(u.des.0)
  j5 <- j4 + 1
  j6 <- ncol(XX) - 1
  
  A.na       <- XX[, j1, drop = FALSE]
  res.na     <- XX[, j2, drop = FALSE]
  u.des.0.na <- XX[, j3:j4, drop = FALSE]
  Z.na       <- XX[, j5:j6, drop = FALSE]
  f.id.na    <- XX[, ncol(XX), drop = FALSE]
  
  active.features <- rowSums(is.na(X)) == 0
  V.na <- V[active.features, active.features]

  # Effective number of observation used for inference (not yet adjusted for df used)
  TT <- nrow(Z.na)
  
  # Remove NA - Out of Sample Uncertainty
  X  <- cbind(e.res, e.des.0)
  XX <- na.omit(X)
  e.res.na   <- XX[, 1, drop = FALSE]
  e.des.0.na <- XX[, -1, drop = FALSE]

  
  ############################################################################
  ############################################################################
  # Proceed cleaning missing data in the post-treatment period
  P.na <- na.omit(P)

  #############################################################################
  ########################################################################
  ## Estimate E[u|H], V[u|H], and Sigma
  # If the model is thought to be misspecified then E[u|H] is estimated
  if (u.missp == TRUE) {
    T.u <- nrow(u.des.0.na) 

    if ((T.u - 20) <= ncol(u.des.0.na)) {
      u.des.0.na <- matrix(1, nrow = T.u, 1)
      if (verbose & (u.order > 0 | u.lags > 0)) {
        warning(paste0("One of u.order > 0 and u.lags > 0 was specified, however the current number of observations (",
                      T.u,") used to estimate conditional moments of the pseudo-residuals ",
                      "is not larger than the number of parameters used in estimation (",ncol(u.des.0.na),") plus 20. ",
                      "To avoid over-fitting issues u.order and u.lags were set to 0."), immediate. = TRUE, call. = FALSE)
      }
      u.order <- 0
      u.lags <- 0
      
    } else {
      
      u.des.0.na <- DUflexGet(u.des.0.na, C, f.id.na, M)
      
    }   
    
    u.mean <- lm(res.na ~ u.des.0.na - 1)$fitted.values 
	params.u <- ncol(u.des.0.na)
    
  } else if (u.missp == FALSE) {
    u.mean <- 0
    params.u <- 0
    T.u <- 0
  }

  # Estimate degrees of freedom to be used for V[u|H] (note that w is pre-regularization)
  df <- df.EST(w.constr = w.constr, w = w, A = A, B = B, J = J, KM = KM)
  
  # Use HC inference to estimate V[u|H]
  result <- u.sigma.est(u.mean = u.mean, u.sigma = u.sigma, res = res.na,
                        Z = Z.na, V = V.na, index = index, TT = TT, M = M, df = df)
  Sigma <- result$Sigma
  Omega <- result$Omega

  Sigma.root <- sqrtm(Sigma)
  
  # Auxiliary logical values to estimate bounds for w
  w.lb.est <- TRUE
  w.ub.est <- TRUE
  if (is.null(w.bounds) == FALSE) {
    if (all(is.na(w.bounds[, 1]) == FALSE)) w.lb.est <- FALSE
    if (all(is.na(w.bounds[, 2]) == FALSE)) w.ub.est <- FALSE
  }
  
  ## Define constrained problem to be simulated
  if (w.lb.est == TRUE | w.ub.est == TRUE) {
    Q <- t(Z.na) %*% V.na %*% Z.na / TT
    colnames(Q) <- colnames(Z.na)
    
    if (w.constr.inf[["p"]] == "no norm") {p.int <- 0}
    if (w.constr.inf[["p"]] == "L1") {p.int <- 1}
    if (w.constr.inf[["p"]] == "L2") {p.int <- 2}
    
    jj <- nrow(P.na)
    
    iters <- round(sims/10)
    perc  <- 0
    # simulate
    if (cores == 1) {
      vsig <- matrix(NA, nrow = sims, ncol = 2*jj)
      
      for (sim in seq_len(sims)) {
        rem <- sim %% iters
        if ((rem == 0) & verbose) {
          perc <- perc + 10
          cat(paste(sim,"/",sims," iterations completed (",perc,"%)"," \r", sep = ""))
          utils::flush.console()
        }
        
        zeta    <- rnorm(length(beta))
        G       <- Sigma.root %*% zeta
        
        for (hor in seq_len(jj)) {
          xt <- P.na[hor,]

          output  <- scpi.in(xt = xt, beta = beta, Q = Q, G = G, J = J, KM = KM, p.int = p.int,
                             QQ = w.constr.inf[["Q"]], lb = w.constr.inf[["lb"]],
                             dire = w.constr.inf[["dir"]], p = w.constr.inf[["p"]], 
                             w.lb.est = w.lb.est, w.ub.est = w.ub.est, opt.list = opt.list.inf)
          
          vsig[sim, hor]      <- output[1]
          vsig[sim, hor + jj] <- output[2]
          
        }
      }
      
    } else if (cores >= 1) {
      
      progress <- function(n) {
        rem <- n %% iters
        if ((rem == 0) & verbose) {
          perc <- n/sims*100
          cat(paste(n,"/",sims," iterations completed (",perc,"%)"," \r", sep = ""))
          utils::flush.console()
        }
      }
      opts <- list(progress=progress)
      
      cl <- parallel::makeCluster(cores)
      doSNOW::registerDoSNOW(cl)
      
      vsig <- foreach::foreach(i = 1 : sims,
                      .packages = c('nloptr','CVXR'),
                      .export   = c('scpi.in','obj.fun.min','obj.fun.max', 'checkConstraints', 'useCVXR',
                                    'single.ineq','double.ineq','norm.equal', 'prepareOptions'),
                      .combine  = rbind,
                      .options.snow = opts) %dorng% {

                      zeta   <- rnorm(length(beta))
                      G      <- Sigma.root %*% zeta
                      
                      ub.sim <- c()
                      lb.sim <- c()
                      
                      for (hor in seq_len(jj)) {
                        xt <- P.na[hor,]
                        
                        output  <- scpi.in(xt = xt, beta = beta, Q = Q, G = G, J = J, KM = KM, p.int = p.int,
                                           QQ = w.constr.inf[["Q"]], lb = w.constr.inf[["lb"]],
                                           dire = w.constr.inf[["dir"]], p = w.constr.inf[["p"]], 
                                           w.lb.est = w.lb.est, w.ub.est = w.ub.est, opt.list = opt.list.inf)

                        lb.sim      <- append(lb.sim, output[1])
                        ub.sim      <- append(ub.sim, output[2])
                        
                      }
                      
                      c(lb.sim, ub.sim)
                    }
      
      parallel::stopCluster(cl)
    }
  }
  
  if (w.lb.est == TRUE) {
    w.lb    <- apply(vsig[, 1:jj, drop = FALSE], 2, quantile, probs = u.alpha/2, na.rm = TRUE)
    fail.lb <- apply(vsig[, 1:jj, drop = FALSE], 2, function(x) sum(is.na(x))/sims*100)
  } else {
    w.lb     <- w.bounds[,1]
    fail.lb  <- rep(0, length(w.bounds[,1]))
  }
    
  if (w.ub.est == TRUE) {
    w.ub    <- apply(vsig[, (jj+1):(2*jj), drop = FALSE], 2, quantile, probs = (1-u.alpha/2), na.rm = TRUE)
    fail.ub <- apply(vsig[, (jj+1):(2*jj), drop = FALSE], 2, function(x) sum(is.na(x))/sims*100)
  } else {
    w.ub     <- w.bounds[,2]
    fail.ub  <- rep(0, length(w.bounds[,2]))
  }
  failed.sims <- rbind(fail.lb, fail.ub)
  rownames(failed.sims) <- c("lb","ub")
  cat("\n")
  
  if ((sum(failed.sims) > 0.1*sims*ncol(vsig)) & verbose) {
    warning("For some of the simulations used to quantify in-sample uncertainty the solution of the optimization problem 
            was not found! We suggest inspecting the magnitude of this issue by consulting the percentage of simulations
            that failed contained in YOUR_SCPI_OBJECT_NAME$inference.results$failed.sims. In case the number of 
            unsuccessful simulations is high, you might want to consider switching the solver or changing the 
            stopping criteria of the algorithm through the option 'opt.list.inf'.", immediate. = TRUE, call. = FALSE)
  }

  ## Adjust for missing values
  P.nomiss <- rowSums(is.na(P)) == 0   # rows of P with no missing values
  
  Wlb <- matrix(NA, nrow = nrow(P), ncol = 1)
  Wub <- matrix(NA, nrow = nrow(P), ncol = 1)
  
  Wlb[P.nomiss, ] <- w.lb
  Wub[P.nomiss, ] <- w.ub
  
  ## PIs for w
  sc.l.0 <- Y.post.fit + Wlb        # Left bound
  sc.r.0 <- Y.post.fit + Wub        # Right bound
  len.0  <- sc.r.0 - sc.l.0         # Length

  
  
  #############################################################################
  #############################################################################
  ## Estimate out-of-sample uncertainty
  #############################################################################
  #############################################################################
  
  # PIs for u
  sc.l.1 <- sc.r.1 <- sc.l.2 <- sc.r.2 <- sc.l.3 <- sc.r.3 <- sc.l.4 <- sc.r.4 <- rep(NA, T1)
  len.1  <- len.2  <- len.3  <- len.4  <- rep(NA, T1)

  # Save mean and variance of e, only for sensitivity analysis
  e.mean <- e.var <- NA

  # Auxiliary logical values to estimate bounds for e
  e.lb.est <- TRUE
  e.ub.est <- TRUE
  if (is.null(e.bounds) == FALSE) {
    if (all(is.na(e.bounds[, 1]) == FALSE)) e.lb.est <- FALSE
    if (all(is.na(e.bounds[, 2]) == FALSE)) e.ub.est <- FALSE
  }
  
  T.e <- nrow(e.des.0.na) 
  
  if ((T.e - 20) <= ncol(e.des.0.na)) {
    e.des.0.na <- matrix(1, nrow = T.e, 1)
    e.des.1 <- matrix(1, nrow = nrow(e.des.1), 1)
    if (verbose & (e.order > 0 | e.lags > 0)) {
      warning(paste0("One of e.order > 0 and e.lags > 0 was specified, however the current number of observations (",
                    T.e,") used to estimate conditional moments of the out-of-sample error",
                    " is not larger than the number of parameters used in estimation (",ncol(e.des.0.na),") plus 20.",
                    " To avoid over-fitting issues e.order and e.lags were set to 0."), immediate. = TRUE, call. = FALSE)
    }
    e.order <- 0
    e.lags <- 0
  }  

  params.e <- ncol(e.des.0.na)

  if (e.method == "gaussian" | e.method == "all") {
    pi.e   <- scpi.out(res = e.res.na, x = e.des.0.na, eval = e.des.1, e.method = "gaussian", alpha = e.alpha/2,
                       e.lb.est = e.lb.est, e.ub.est =  e.ub.est)
    
    # Overwrite with user's input
    if (e.lb.est == FALSE) pi.e$lb <- e.bounds[, 1]
    if (e.ub.est == FALSE) pi.e$ub <- e.bounds[, 2]
    
    sc.l.1 <- sc.l.0 + pi.e$lb
    sc.r.1 <- sc.r.0 + pi.e$ub
    len.1  <- sc.r.1 - sc.l.1
    e.mean <- pi.e$e.1
    e.var  <- pi.e$e.2
  }

  if (e.method == "ls" | e.method == "all") {
    pi.e   <- scpi.out(res = e.res.na, x = e.des.0.na, eval = e.des.1, e.method = "ls", alpha = e.alpha/2,
                       e.lb.est = e.lb.est, e.ub.est =  e.ub.est)
    
    # Overwrite with user's input
    if (e.lb.est == FALSE) pi.e$lb <- e.bounds[, 1]
    if (e.ub.est == FALSE) pi.e$ub <- e.bounds[, 2]

    sc.l.2 <- sc.l.0 + pi.e$lb
    sc.r.2 <- sc.r.0 + pi.e$ub
    len.2  <- sc.r.2 - sc.l.2
    e.mean <- pi.e$e.1
    e.var  <- pi.e$e.2
  }

  if (e.method == "qreg" | e.method == "all") {
    if (e.order == 0) {
      
      lb <- quantile(e.res.na, e.alpha/2)
      ub <- quantile(e.res.na, 1-e.alpha/2)
      
      # Overwrite with user's input
      if (e.lb.est == FALSE) lb <- e.bounds[, 1]
      if (e.ub.est == FALSE) ub <- e.bounds[, 2]
      
      sc.l.3 <- sc.l.0 + lb
      sc.r.3 <- sc.r.0 + ub
      len.3  <- sc.r.3 - sc.l.3
      
    } else {
      pi.e   <- scpi.out(res = e.res.na, x = e.des.0.na, eval = e.des.1, e.method = "qreg", alpha = e.alpha/2,
                         e.lb.est = e.lb.est, e.ub.est =  e.ub.est)

      # Overwrite with user's input
      if (e.lb.est == FALSE) pi.e$lb <- e.bounds[, 1]
      if (e.ub.est == FALSE) pi.e$ub <- e.bounds[, 2]
      
      sc.l.3 <- sc.l.0 + pi.e$lb
      sc.r.3 <- sc.r.0 + pi.e$ub
      len.3  <- sc.r.3 - sc.l.3
    }
  }


  #############################################################################
  #############################################################################
  ## Return objects

  CI.0 <- cbind(sc.l.0, sc.r.0, len.0)
  colnames(CI.0) <- c("Left Bound", "Right Bound", "Length")
  

  if (is.null(sc.l.1) == FALSE) {
    CI.1 <- cbind(sc.l.1, sc.r.1, len.1)
    colnames(CI.1) <- c("Left Bound", "Right Bound", "Length")
  } else {
    CI.1 <- NULL
  }

  if (is.null(sc.l.2) == FALSE) {
    CI.2 <- cbind(sc.l.2, sc.r.2, len.2)
    colnames(CI.2) <- c("Left Bound", "Right Bound", "Length")
  } else {
    CI.2 <- NULL
  }

  if (is.null(sc.l.3) == FALSE) {
    CI.3 <- cbind(sc.l.3, sc.r.3, len.3)
    colnames(CI.3) <- c("Left Bound", "Right Bound", "Length")
  } else {
    CI.3 <- NULL
  }
  
  u.user <- is.null(u.design) == FALSE # True if user provided the design matrix for in-sample inference
  e.user <- is.null(e.design) == FALSE # True if user provided the design matrix for out-of-sample inference

  inference.results <- list(  CI.in.sample    = CI.0,
                              CI.all.gaussian = CI.1,
                              CI.all.ls       = CI.2,
                              CI.all.qreg     = CI.3,
                              Sigma           = Sigma,
                              u.mean          = u.mean,
                              u.var           = Omega,
                              e.mean          = e.mean,
                              e.var           = e.var,
                              u.missp         = u.missp,
                              u.lags          = u.lags,
                              u.order         = u.order,
                              u.sigma         = u.sigma,
                              u.user          = u.user,
                              u.T             = T.u,
                              u.params        = params.u,
                              u.D             = u.des.0.na,
                              e.method        = e.method,
                              e.lags          = e.lags,
                              e.order         = e.order,
                              e.user          = e.user,
                              e.T             = T.e,
                              e.params        = params.e,
                              e.D             = e.des.0.na,
                              rho             = rho,
                              u.alpha         = u.alpha,
                              e.alpha         = e.alpha,
                              sims            = sims,
                              failed.sims     = failed.sims)


  result <- list( data               = sc.pred$data,
                  est.results        = sc.pred$est.results,
                  inference.results  = inference.results)

  class(result) <- 'scpi_scpi'

  #############################################################################
  #############################################################################
  ## Plot
  
  if (plot == TRUE) {
    if (is.null(plot.name) == FALSE) {
      fig.name = plot.name
    } else {
      fig.name = "scpi_default_plot"
    }
    
     scplot(result = result, fig.path = getwd(),
            fig.name = fig.name, fig.format = "png", save.data = save.data)
  }



    return(result)
}

