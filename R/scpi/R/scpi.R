###############################################################################

#' @title Prediction Intervals for Synthetic Control Methods
#'
#' @description The command implements estimation and inference procedures for Synthetic Control (SC) methods using least squares, lasso, ridge, or simplex-type constraints. Uncertainty is quantified using prediction
#' intervals according to \href{https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, Feng, and Titiunik (2021)}. \code{\link{scpi}} returns the estimated 
#' post-treatment series for the synthetic unit through the command \code{\link{scest}} and quantifies in-sample and out-of-sample uncertainty to provide confidence intervals
#' for each point estimate.
#'
#' Companion \href{https://www.stata.com/}{Stata} and \href{https://www.python.org/}{Python} packages are described in \href{https://arxiv.org/abs/2202.05984}{Cattaneo, Feng, Palomba, and Titiunik (2022)}.
#'
#' Companion commands are:  \link{scdata} and \link{scdataMulti} for data preparation in the single and multiple treated unit(s) cases, respectively,
#' \link{scest} for point estimation, \link{scplot} and \link{scplotMulti} for plots in the single and multiple treated unit(s) cases, respectively.
#'
#' Related Stata, R, and Python packages useful for inference in SC designs are described in the following website:
#'
#' \href{ https://nppackages.github.io/scpi/}{ https://nppackages.github.io/scpi/}
#'
#' For an introduction to synthetic control methods, see \href{https://www.aeaweb.org/articles?id=10.1257/jel.20191450}{Abadie (2021)} and references therein.
#'
#' @param data a class 'scdata' object, obtained by calling \code{\link{scdata}}, or class 'scdataMulti' obtained via \code{\link{scdataMulti}}.
#' @param w.constr a list specifying the constraint set the estimated weights of the donors must belong to.
#' \code{w.constr} can contain up to five elements:
#' - `\code{p}', a scalar indicating the norm to be used (\code{p} should be one of "no norm", "L1", and "L2")
#' - `\code{dir}', a string indicating whether the constraint on the norm is an equality ("==") or inequality ("<=")
#' - `\code{Q}', a scalar defining the value of the constraint on the norm
#' - `\code{lb}', a scalar defining the lower bound on the weights. It can be either 0 or \code{-Inf}.
#' - `\code{name}', a character selecting one of the default proposals
#' See the \strong{Details} section for more.
#' @param P a \eqn{I\cdot T_1\times I\cdot (J+KM)} matrix containing the design matrix to be used to obtain the predicted.
#' post-intervention outcome of the synthetic control unit. \eqn{T_1} is the number of post-treatment periods,
#' \eqn{J} is the size of the donor pool, and \eqn{K_1} is the number of covariates used for adjustment in the outcome equation.
#' @param V specifies the weighting matrix to be used when minimizing the sum of squared residuals
#' \deqn{(\mathbf{A}-\mathbf{B}\mathbf{w}-\mathbf{C}\mathbf{r})'\mathbf{V}(\mathbf{A}-\mathbf{B}\mathbf{w}-\mathbf{C}\mathbf{r})}
#' The default is the identity matrix, so equal weight is given to all observations. In the case of multiple treated observations
#' (you used \code{\link{scdataMulti}} to prepare the data), the user can specify \code{V} as a string equal to either "separate" or "pooled".
#' In both cases, the user can provide a conformable matrix as input. See the \strong{Details} section for more.
#' @param rho a string specifying the regularizing parameter that imposes sparsity on the estimated vector of weights. If
#' \code{rho = 'type-1'} (the default), then the tuning parameter is computed based on optimization inequalities. Users can provide a scalar 
#' with their own value for \code{rho}. Other options are described in the \strong{Details} section.
#' @param rho.max a scalar indicating the maximum value attainable by the tuning parameter \code{rho}.
#' @param lgapp selects the way local geometry is approximated in simulation. The options are "generalized"
#' and "linear". The first one accommodates for possibly non-linear constraints, whilst the second one is valid
#' with linear constraints only.
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
#' @param plot a logical specifying whether \code{\link{scplot}} should be called and a plot saved in the current working
#' directory. For more options see \code{\link{scplot}}.
#' @param plot.name a string containing the name of the plot (the format is by default .png). For more options see \code{\link{scplot}}.
#' @param cores number of cores to be used by the command. The default is one.
#' @param w.bounds a \eqn{N_1\cdot T_1\times 2} matrix with the user-provided bounds on \eqn{\beta}. If \code{w.bounds} is provided, then
#' the quantification of in-sample uncertainty is skipped. It is possible to provide only the lower bound or the upper bound
#' by filling the other column with \code{NA}s.
#' @param e.bounds a \eqn{N_1\cdot T_1\times 2} matrix with the user-provided bounds on \eqn{(\widehat{\mathbf{w}},
#' \widehat{\mathbf{r}})^{\prime}}. If \code{e.bounds} is provided, then
#' the quantification of out-of-sample uncertainty is skipped. It is possible to provide only the lower bound or the upper bound
#' by filling the other column with \code{NA}s.
#'
#' @param save.data a character specifying the name and the path of the saved dataframe containing the processed data used to produce the plot. 
#'
#' @param verbose if \code{TRUE} prints additional information in the console.
#'
#' @return
#' The function returns an object of class 'scpi' containing three lists. The first list is labeled 'data' and contains used
#' data as returned by \code{\link{scdata}} and some other values.
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
#' \item{\code{donors.list}, list containing the IDs of the donors of each treated unit}}}
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
#' \item{bounds}{a list containing the estimated bounds (in-sample and out-of-sample uncertainty).}
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
#' \item{Q.star}{a list containing the regularized constraint on the norm.}
#' \item{epskappa}{a vector containing the estimates for \eqn{\epsilon_{\kappa}}.}
#' \item{sims}{an integer indicating the number of simulations used in quantifying in-sample uncertainty.}
#' \item{failed.sims}{a matrix containing the number of failed simulations per post-treatment period to estimate lower and upper bounds.}
#'
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
#' where \eqn{Q} is a tuning parameter that is by default computed as
#' \deqn{(J+KM) \widehat{\sigma}_u^{2}/||\widehat{\mathbf{w}}_{OLS}||_{2}^{2}}
#' where \eqn{J} is the number of donors and \eqn{KM} is the total number of covariates used for adjustment.
#' The user can provide \code{Q} as an element of the list (eg. \code{w.constr =
#' list(name = "ridge", Q = 1)}).}
#'
#' \item{If \code{name == "ols"}, then the problem is unconstrained and the vector of weights
#' is estimated via ordinary least squares.}
#' 
#' \item{If \code{name == "L1-L2"}, then
#' \deqn{||\mathbf{w}||_1 = 1,\:\:\: ||\mathbf{w}||_2 \leq Q,}
#' where \eqn{Q} is a tuning parameter computed as in the "ridge" case.}
#' }}
#'
#' \item{\strong{Weighting Matrix.}
#' \itemize{
#' \item{if \code{V <- "separate"}, then \eqn{\mathbf{V} = \mathbf{I}} and the minimized objective function is
#' \deqn{\sum_{i=1}^{N_1} \sum_{l=1}^{M} \sum_{t=1}^{T_{0}}\left(a_{t, l}^{i}-\mathbf{b}_{t, l}^{{i \prime }} \mathbf{w}^{i}-\mathbf{c}_{t, l}^{{i \prime}} \mathbf{r}_{l}^{i}\right)^{2},}
#' which optimizes the separate fit for each treated unit.}
#' \item{if \code{V <- "pooled"}, then \eqn{\mathbf{V} = \mathbf{1}\mathbf{1}'\otimes \mathbf{I}} and the minimized objective function is
#' \deqn{\sum_{l=1}^{M} \sum_{t=1}^{T_{0}}\left(\frac{1}{N_1^2} \sum_{i=1}^{N_1}\left(a_{t, l}^{i}-\mathbf{b}_{t, l}^{i \prime} \mathbf{w}^{i}-\mathbf{c}_{t, l}^{i\prime} \mathbf{r}_{l}^{i}\right)\right)^{2},}
#' which optimizes the pooled fit for the average of the treated units.}
#' \item{if \code{V} is a user-provided matrix, then in must be a \eqn{v\times v} positive-definite matrix where \eqn{v} is the 
#' number of rows of \eqn{\mathbf{B}} (or \eqn{\mathbf{C}}) after missing values have been taken into account. In case the user
#' wants to provide their own \code{V}, we suggest to check the appropriate dimension \eqn{v} by inspecting the output
#' of either \code{scdata} or \code{scdataMulti}.}
#' }}
#'
#'
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
#' then \eqn{\mathbf{E}[\mathbf{u}\: |\: \mathbf{D}_u]} is estimated using a linear regression of
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
#' number of observation used to estimate moments of \eqn{\mathbf{u}} is \eqn{N_1\cdot T_0\cdot M}, whilst for moments of \eqn{e_t} is \eqn{T_0}.
#' Our rule of thumb to avoid over-fitting is to check if \eqn{N_1\cdot T_0\cdot M \geq d_u + 10} or \eqn{T_0 \geq d_e + 10}. If the 
#' former condition is not satisfied we automatically set \code{u.order = u.lags = 0}, if instead the latter is not met
#' we automatically set \code{e.order = e.lags = 0}.}
#'
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
#' \item{\href{https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, M. D., Feng, Y., and Titiunik, R. 
#' (2021)}. Prediction intervals for synthetic control methods. \emph{Journal of the American Statistical Association}, 116(536), 1865-1880.}
#' \item{\href{https://arxiv.org/abs/2202.05984}{Cattaneo, M. D., Feng, Y., Palomba F., and Titiunik, R. (2022).}
#' scpi: Uncertainty Quantification for Synthetic Control Methods, \emph{arXiv}:2202.05984.}
#' \item{\href{https://arxiv.org/abs/2210.05026}{Cattaneo, M. D., Feng, Y., Palomba F., and Titiunik, R. (2022).}
#' Uncertainty Quantification in Synthetic Controls with Staggered Treatment Adoption, \emph{arXiv}:2210.05026.}
#' }
#'
#' @seealso \code{\link{scdata}}, \code{\link{scdataMulti}}, \code{\link{scest}}, \code{\link{scplot}}, \code{\link{scplotMulti}}
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
#' result <- scpi(df, w.constr = list(name = "simplex", Q = 1), cores = 1, sims = 10)
#' result <- scpi(df, w.constr = list(lb = 0, dir = "==", p = "L1", Q = 1),
#'                cores = 1, sims = 10)
#'                            
#' @export

scpi  <- function(data,
                  w.constr     = NULL,
                  V            = "separate",
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
                  rho.max      = 0.2,
                  lgapp        = "generalized",
                  cores        = 1,
                  plot         = FALSE,
                  plot.name    = NULL,
                  w.bounds     = NULL,
                  e.bounds     = NULL,
                  save.data    = NULL,
                  verbose      = TRUE) {


  if ( (methods::is(data, "scdata") || methods::is(data, "scdataMulti")) == FALSE ) {
    stop("data should be the object returned by running scdata or scdata_multi!")
  }

  if (methods::is(data, 'scdata') == TRUE) {
    class.type <- 'scpi_data'
  } else if (methods::is(data, 'scdataMulti') == TRUE) {
    class.type <- 'scpi_data_multi'
  }

  V.type <- V

  #############################################################################
  #############################################################################
  ## Estimation of synthetic weights
  if (verbose) {
    cat("---------------------------------------------------------------\n")
    cat("Estimating Weights...\n")
  }
  sc.pred <- scest(data = data, w.constr = w.constr, V = V)


  #############################################################################
  #############################################################################
  ## Retrieve processed data from scest

  A           <- sc.pred$data$A                           # Features of treated unit
  B           <- sc.pred$data$B                           # Features of control units
  C           <- sc.pred$data$C                           # Covariates for adjustment
  Z           <- sc.pred$data$Z                           # B and C column-bind
  Y.donors    <- data$Y.donors                            # Outcome variable of control units
  K           <- sc.pred$data$specs$K                     # Number of covs for adjustment per feature
  KM          <- sc.pred$data$specs$KM                    # Dimension of r (total number of covs for adj)
  J           <- sc.pred$data$specs$J                     # Number of donors
  M           <- sc.pred$data$specs$M                     # Number of features
  T0          <- sc.pred$data$specs$T0.features           # Time periods used per feature
  T1          <- sc.pred$data$specs$T1.outcome            # Number of out-of-sample periods
  features    <- sc.pred$data$specs$features              # Name of features
  constant    <- sc.pred$data$specs$constant              # Logical indicating whether a constant is included
  out.feat    <- sc.pred$data$specs$out.in.features       # Logical indicating whether the outcome variable is among features
  coig.data   <- sc.pred$data$specs$cointegrated.data     # Logical indicating whether B is cointegrated
  w.constr    <- sc.pred$est.results$w.constr             # Constraints on w
  V           <- sc.pred$est.results$V                    # Weighting matrix
  w           <- sc.pred$est.results$w                    # Estimated vector of weights
  r           <- sc.pred$est.results$r                    # Estimated coefficients of covariates
  b           <- sc.pred$est.results$b                    # w and r column-bind
  Y.post.fit  <- sc.pred$est.results$Y.post.fit           # Estimated post-treatment outcome for SC unit
  res         <- sc.pred$est.results$res                  # Residuals from estimation
  outcome.var <- sc.pred$data$specs$outcome.var           # name of outcome variable
  sc.effect   <- sc.pred$data$specs$effect                # Causal quantity of interest
  sparse.mat  <- sc.pred$data$specs$sparse.matrices       # Whether sparse matrices are involved or not
  
  if (class.type == 'scpi_data') {
    Jtot            <- J
    KMI             <- KM
    I               <- 1
    T0.tot          <- sum(T0)                            # Total number of observations used in estimation
    T0.M            <- T0.tot
    T1.tot          <- T1
    features        <- list(features)
    out.feat        <- list(out.feat)
    T0              <- list(T0)
    names(T0)       <- sc.pred$data$specs$treated.units

  } else if (class.type == 'scpi_data_multi') {
    J               <- unlist(J)
    Jtot            <- sum(J)
    KMI             <- data$specs$KMI                              # total number of covariates used for adjustment
    I               <- data$specs$I                                # number of treated units
    T0.M            <- unlist(lapply(data$specs$T0.features, sum)) # observations per treated unit
    T0.tot          <- sum(T0.M)                                   # Total number of observations used in estimation
    T1.tot          <- sum(unlist(T1))                             # Total number of observations post-treatment
  }

  if (sc.effect == "unit") {
    T1.tot <- I
    T1 <- lapply(T1, function(x) 1)
  }

  # Check on P
  if (is.null(P) == TRUE) {
    P <- sc.pred$data$P                         # Matrix for out-of-sample prediction

  } else { # User-provided prediction matrix P (should be T1 by (J+KM))
    if (is.matrix(P) == FALSE) {
      stop("The object P should be a matrix!")
    }
    if (class.type == 'scpi_data') {
      Prows <- T1
      Pcols <- J+KM
    } else if (class.type == 'scpi_data_multi') {
      Prows <- T1.tot
      Pcols <- Jtot + KMI
    }

    if (nrow(P) != Prows) {
      stop(paste("The matrix P currently has", nrow(P), "rows when instead", Prows, "were expected
                 (i.e. the number of post-intervention periods)!"))
    }
    if (ncol(P) != Pcols) {
      stop(paste("The matrix P currently has", ncol(P), "columns when instead", Pcols, "were expected
                 (i.e. the size of the donor pool plus the number of covariates used in adjustment in the outcome equation)!"))
    }
  }

  if (!(e.method %in% c("gaussian", "ls", "qreg", "all"))) {
    stop("The object e.method should be one of 'gaussian', 'ls', 'qreg', or 'all'.")
  }

  if (!(lgapp %in% c("linear", "generalized"))) {
    stop("The object lgapp should be one of 'linear' or 'generalized'.")
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
      stop(paste0("w.bounds should be a matrix with ", length(Y.post.fit),
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
      if (!(rho %in% c('type-1','type-2','type-3'))) {
        stop("When not a scalar, 'rho' must be 'type-1', 'type-2', or 'type-3'.")
      }
    }
  } else {
    rho <- "type-1"
  }

  # Check on number of cores
  if (is.null(cores) == FALSE) {
    n.cores <- parallel::detectCores(logical = TRUE)
    if (cores > n.cores) {
      stop(paste("You selected", cores, "cores, but only", n.cores, " cores are available on your machine!"))
    }
  } else {
    cores <- parallel::detectCores(logical = TRUE) - 1
    warning(paste("scpi is using",cores,"cores for estimation! You can adjust the number of cores with the 'cores' option."), immediate. = TRUE, call. = FALSE)
  }

  if (verbose) {
    if (class.type == 'scpi_data') {
      constr.type <- w.constr[['name']]
    } else if (class.type == 'scpi_data_multi') {
      constr.type <- w.constr[[1]]$name
    }

    cat("Quantifying Uncertainty\n")
    #executionTime(T0.tot, Jtot, I, T1.tot, sims, cores, constr.type)
  }

  # create lists of matrices
  A.list   <- mat2list(A)
  B.list   <- mat2list(B)
  if (is.null(C) == FALSE) {
    if (ncol(C) > 0) {
      C.list   <- mat2list(C)
    } else {
      C.list <- rep(list(NULL), I)
      names(C.list) <- sc.pred$data$specs$treated.units
    }
  } else {
    C.list <- rep(list(NULL), I)
    names(C.list) <- sc.pred$data$specs$treated.units
  }

  if (sc.pred$data$specs$effect == "time") {
    P.list <- list()
    names <- strsplit(colnames(P), "\\.")
    cnames <- unlist(lapply(names, "[[", 1))
    for (tr in sc.pred$data$specs$treated.units) {
      P.list[[tr]] <- P[, cnames == tr, drop=FALSE]
    }
  } else {
    P.list <- mat2list(P)
  }

  if (!is.null(sc.pred$data$P.diff)) {
    Pd.list <- mat2list(sc.pred$data$P.diff)
  } else {
    Pd.list <- rep(list(NULL), I)
  }
  
  V.list   <- mat2list(V)
  w.list   <- mat2list(as.matrix(w))
  res.list <- mat2list(res)
  Y.d.list <- mat2list(Y.donors)

  if (sparse.mat == TRUE) {
    res.list <- lapply(res.list, as.matrix)
    B.list <- lapply(B.list, as.matrix)
    if (is.null(C) == FALSE) C.list <- lapply(C.list, as.matrix)
    P.list <- lapply(P.list, as.matrix)
  }
  
  #############################################################################
  #############################################################################
  ### Estimate In-Sample Uncertainty
  #############################################################################
  #############################################################################

  if (class.type == 'scpi_data') {
    w.constr.list <- list(w.constr)
    names(w.constr.list) <- sc.pred$data$specs$treated.units
  } else if (class.type == 'scpi_data_multi') {
    w.constr.list <- w.constr
  }

  w.star <- index.w <- rho.vec <- Q.star <- Q2.star <- f.id <- e.res <- u.names <- e.rownames <- e.colnames <- e1.rownames <- c()
  u.des.0 <- e.des.0 <- e.des.1 <- matrix(NA, 0, 0)
  w.constr.inf <- list()

  for (i in seq_len(I)) {
    ## Regularize W and local geometry (treated unit by treated unit)
    loc.geom <- local.geom(w.constr.list[[i]], rho, rho.max, res.list[[i]], B.list[[i]], 
                           C.list[[i]], coig.data[[i]], T0.M[[i]], J[[i]], w.list[[i]],
                           verbose)
    
    w.star       <- c(w.star, loc.geom$w.star)
    index.w      <- c(index.w, loc.geom$index.w)
    w.constr.inf <- append(w.constr.inf, list(loc.geom$w.constr))
    rho.vec      <- c(rho.vec, loc.geom$rho)
    Q.star       <- c(Q.star, loc.geom$Q.star)
    Q2.star      <- c(Q2.star, loc.geom$Q2.star)
    index.i      <- c(loc.geom$index.w, rep(TRUE, KM[[i]]))
    
    # Extract feature id from rownames of B
    feature.id <- unlist(purrr::map(stringr::str_split(rownames(B.list[[i]]), "\\."), 2))

    for (f in features[[i]]) {
      len.feat <- sum(feature.id == f)
      if (len.feat <= u.lags && u.lags > 0) {
        if (verbose) {
          warning("At least one of your features is observed for less periods than the number of lags, u.lags reverted to 0.", immediate. = TRUE, call. = FALSE)
        }
        u.lags <- 0
       }
    }

    ## Prepare design matrix for in-sample uncertainty
    obj <- u.des.prep(B.list[[i]], C.list[[i]], u.order, u.lags, coig.data[[i]],
                      T0.M[i], constant[[i]], index.i, loc.geom$index.w,
                      features[[i]], feature.id, u.design, res.list[[i]])
    u.names <- c(u.names, colnames(obj$u.des.0))
    
    u.des.0 <- Matrix::bdiag(u.des.0, obj$u.des.0)
    f.id <- c(f.id, as.factor(feature.id))
    
    ## Prepare design matrices for out-of-sample uncertainty
    e.des <- e.des.prep(B.list[[i]], C.list[[i]], P.list[[i]], e.order, e.lags,
                        res.list[[i]], sc.pred, Y.d.list[[i]], out.feat[[i]],
                        features[[i]], J[[i]], index.i, loc.geom$index.w,
                        coig.data[[i]], T0[[i]][outcome.var], T1[[i]], constant[[i]], 
                        e.design, Pd.list[[i]], sc.pred$data$specs$effect, I, class.type)
 
    e.res   <- c(e.res, e.des$e.res)
    e.rownames <- c(e.rownames, rownames(e.des$e.res))
    cnames <- rep(paste0(names(w.constr.list)[[i]], "."), ncol(e.des$e.des.0))
    e.colnames <- c(e.colnames, cnames)

    if (sc.pred$data$specs$effect == "time") {
      trname <- unlist(purrr::map(stringr::str_split(rownames(e.des$e.des.0)[1], "\\."), 1))
      rnames <- paste(trname, as.character(c(1:nrow(e.des$e.des.1))), sep=".")
      e1.rownames <- c(e1.rownames, rnames)
    }
    
    e.des.0 <- Matrix::bdiag(e.des.0, e.des$e.des.0)
    e.des.1 <- Matrix::bdiag(e.des.1, e.des$e.des.1)
  }

  # Create an index that selects all non-zero weights and additional covariates
  index <- c(index.w, rep(TRUE, KMI))

  if (lgapp == "generalized") {
    beta <- b # we use rho only to impose sparsity on B when predicting moments
    Q <- c()
    for (i in seq_len(I)) {
      Q <- c(Q, w.constr.list[[i]]$Q)
    }

    reg.geom <- local.geom.2step(w, r, rho.vec, w.constr.list, Q, I)
    Q.star <- reg.geom$Q
    lb <- reg.geom$lb
    
  } else if (lgapp == "linear") { # we use rho to regularize w too
    beta  <- c(w.star, r)
    lb <- c()
    for (i in seq_len(I)) {
      lb <- c(lb, rep(w.constr.inf[[i]]$lb, J[i]))
    }
  }

  names(beta) <- names(sc.pred$est.results$b)

  # Transform sparse matrices to matrices
  e.des.0 <- as.matrix(e.des.0)
  e.des.1 <- as.matrix(e.des.1)
  u.des.0 <- as.matrix(u.des.0)
  e.res   <- as.matrix(e.res)
  colnames(u.des.0) <- u.names
  rownames(e.res) <- e.rownames
  rownames(e.des.0) <- e.rownames
  colnames(e.des.0) <- e.colnames
  if (sc.pred$data$specs$effect == "time")  {
    rownames(e.des.1) <- e1.rownames
  } else {
    rownames(e.des.1) <- rownames(P)
  }
  colnames(e.des.1) <- e.colnames

  #############################################################################
  ###########################################################################
  # Remove NA - In Sample Uncertainty
  X  <- cbind(A, res, u.des.0, Z, f.id)
  # na.omit does not work with sparse matrices, so have to do it manually
  rowKeep <- Matrix::rowSums(is.na(X)) == 0
  XX <- X[rowKeep, ]
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

  active.features <- Matrix::rowSums(is.na(X)) == 0
  V.na <- V[active.features, active.features]

  # Effective number of observation used for inference (not yet adjusted for df used)
  TT <- nrow(Z.na)

  # Remove NA - Out of Sample Uncertainty
  X  <- cbind(e.res, e.des.0)
  rowKeep <- Matrix::rowSums(is.na(X)) == 0
  XX <- X[rowKeep, ]
  e.res.na   <- XX[, 1, drop = FALSE]
  e.des.0.na <- XX[, -1, drop = FALSE]

  # Proceed cleaning missing data in the post-treatment period
  rowKeep <- Matrix::rowSums(is.na(P)) == 0
  P.na <- P[rowKeep, ]

  #############################################################################
  ########################################################################
  ## Estimate E[u|H], V[u|H], and Sigma
  # If the model is thought to be misspecified then E[u|H] is estimated
  
  if (u.missp == TRUE) {
    T.u <- nrow(u.des.0.na)
    u.des.list <- mat2list(as.matrix(u.des.0.na)) # as.matrix takes care of the sparse case
    f.id.list <- mat2list(as.matrix(f.id.na))
    u.des.0.flex <- matrix(NA, 0, 0)
    u.des.0.noflex <- matrix(NA, 0, 0)
    for (i in seq_len(I)) {
      u.des.0.flex <- Matrix::bdiag(u.des.0.flex, 
                                    DUflexGet(as.matrix(u.des.list[[i]]), C.list[[i]], f.id.list[[i]], M[[i]]))
      u.des.0.noflex <- Matrix::bdiag(u.des.0.noflex, as.matrix(u.des.list[[i]]))
    }

    df.U <- T.u - 10

    u.simple <- df.U <= ncol(u.des.0.noflex)
    u.noflex <- (df.U > ncol(u.des.0.noflex)) & (df.U <= ncol(u.des.0.flex))
    u.flex <- df.U > ncol(u.des.0.flex)

    if (u.simple) {
      if (verbose && (u.order > 0 || u.lags > 0)) {
        warning(paste0("One of u.order > 0 and u.lags > 0 was specified, however the current number of observations (",
                      T.u, ") used to estimate conditional moments of the pseudo-residuals ",
                      "is not larger than the number of parameters used in estimation (",ncol(u.des.0.flex),") plus 10. ",
                      "To avoid over-fitting issues u.order and u.lags were set to 0."), immediate. = TRUE, call. = FALSE)
      }
      u.des.0.na <- matrix(1, nrow = T.u, 1)
      u.order <- 0
      u.lags <- 0
    } else if (u.noflex) {
      if (verbose && (u.order > 0 || u.lags > 0)) {
        warning(paste0("The current number of observations (",T.u,") used to estimate conditional moments of the pseudo-residuals ",
                       "is not larger than the number of parameters used in estimation (",ncol(u.des.0.flex),") plus 10 when allowing for a ",
                       "feature specific model. To avoid over-fitting issues, the conditional moments of the pseudo-residuals are predicted with ",
                       "the same model across features."), immediate. = TRUE, call. = FALSE)
      }
      u.des.0.na <- as.matrix(u.des.0.noflex)

    } else if (u.flex){
      u.des.0.na <- as.matrix(u.des.0.flex)

    }

    y <- res.na[, 1] # takes care of the sparse matrices
    u.mean <- lm(y ~ u.des.0.na - 1)$fitted.values
    params.u <- ncol(u.des.0.na)

  } else if (u.missp == FALSE) {
    u.mean <- 0
    params.u <- 0
    T.u <- 0
  }

  # Estimate degrees of freedom to be used for V[u|H] (note that w is pre-regularization)
  df <- df.EST(w.constr = w.constr.list[[1]], w = w, B = B, J = Jtot, KM = KMI)
  if (df >= TT) {
    df <- TT - 1
    warning(paste0("The current specification uses more degrees of freedom than observations. We suggest to increase the level of ",
                   "sparsity or consider using a smaller donor pool."), immediate. = TRUE, call. = FALSE)
  }

  # Use HC inference to estimate V[u|H]
  result <- u.sigma.est(u.mean = u.mean, u.sigma = u.sigma, res = as.matrix(res.na),
                        Z = Z.na, V = V.na, index = index, TT = TT, df = df)

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
  if (w.lb.est == TRUE || w.ub.est == TRUE) {
    vsigg <- insampleUncertaintyGet(Z.na, V.na, P.na, beta, Sigma.root, J, KMI, I,
                                    w.constr.inf[[1]], Q.star, Q2.star, lb, TT, sims, cores, verbose, 
                                    w.lb.est, w.ub.est)

    vsig <- vsigg[rowSums(is.na(vsigg)) < ncol(vsigg), ] # remove simulations were SOCP was not solved at any horizon
  }

  if (w.lb.est == TRUE) {
    w.lb    <- apply(vsig[, 1:nrow(P.na), drop = FALSE], 2, quantile, probs = u.alpha / 2, na.rm = TRUE)
    fail.lb <- apply(vsigg[, 1:nrow(P.na), drop = FALSE], 2, function(x) sum(is.na(x)) / sims * 100)

  } else {
    w.lb     <- w.bounds[, 1]
    fail.lb  <- rep(0, length(w.bounds[, 1]))

  }

  if (w.ub.est == TRUE) {
    w.ub    <- apply(vsig[, (nrow(P.na) + 1):(2 * nrow(P.na)), drop = FALSE], 2, quantile,
                     probs = (1 - u.alpha / 2), na.rm = TRUE)
    fail.ub <- apply(vsigg[, (nrow(P.na) + 1):(2 * nrow(P.na)), drop = FALSE], 2, function(x) sum(is.na(x))/sims*100)

  } else {
    w.ub     <- w.bounds[, 2]
    fail.ub  <- rep(0, length(w.bounds[, 2]))

  }

  if (w.lb.est == FALSE && w.ub.est == FALSE) {
    vsig <- matrix(0, nrow = sims, ncol = 2 * T1.tot)
  }
  failed.sims <- rbind(fail.lb, fail.ub)
  rownames(failed.sims) <- c("lb", "ub")
  cat("\n")

  if ((sum(failed.sims) > 0.1 * sims * ncol(vsig)) && verbose) {
    warning("For some of the simulations used to quantify in-sample uncertainty the solution of the optimization problem 
          was not found! We suggest inspecting the magnitude of this issue by consulting the percentage of simulations
          that failed contained in YOUR_SCPI_OBJECT_NAME$inference.results$failed.sims.",
            immediate. = TRUE, call. = FALSE)
  }

  ## Adjust for missing values
  P.nomiss <- Matrix::rowSums(is.na(P)) == 0   # rows of P with no missing values

  Wlb <- matrix(NA, nrow = nrow(P), ncol = 1)
  Wub <- matrix(NA, nrow = nrow(P), ncol = 1)
  rownames(Wlb) <- rownames(P)
  rownames(Wub) <- rownames(P)

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
  if (w.constr.inf[[1]][["p"]] %in% c("L2","L1-L2")) {
    beta.mat <- as.matrix(beta)
    epsk <- epskappaGet(P, rho.vec, beta.mat, I, sc.pred$data$specs$effect)
    epsk.j <- epskappaGet(P, rho.vec, beta.mat, I, sc.pred$data$specs$effect, joint = TRUE)
  } else {
    epsk <- epsk.j <- 0
  }

  # PIs for u
  sc.l.1 <- sc.r.1 <- sc.l.2 <- sc.r.2 <- sc.l.3 <- sc.r.3 <- sc.l.4 <- sc.r.4 <- rep(NA, T1.tot)
  len.1  <- len.2  <- len.3  <- len.4  <- rep(NA, T1.tot)

  # Save mean and variance of e, only for sensitivity analysis
  e.mean <- e.var <- NA

  # Auxiliary logical values to estimate bounds for e
  e.lb.est <- TRUE
  e.ub.est <- TRUE
  if (is.null(e.bounds) == FALSE) {
    if (all(is.na(e.bounds[, 1]) == FALSE)) e.lb.est <- FALSE
    if (all(is.na(e.bounds[, 2]) == FALSE)) e.ub.est <- FALSE
  }

  if (sc.pred$data$specs$effect == "time") {
    scale.x <- sc.pred$data$specs$I
  } else {
    scale.x <- 1
  }

  T.e <- nrow(e.des.0.na)
  params.e <- ncol(e.des.0.na)
  if ((T.e - 10) <= params.e) {
    names0 <- rownames(e.des.0.na)
    names1 <- rownames(e.des.1)

    e.des.0.na.list <- mat2list(e.des.0.na)
    e.des.1.list <- mat2list(e.des.1)

    for (i in seq_len(I)) {
      if (sc.effect == "time") {
        aux <- matrix(1/scale.x, nrow = nrow(e.des.0.na.list[[i]]), 1)
        auxx <- matrix(1/scale.x, nrow = nrow(e.des.1.list[[i]]), 1)
      } else {
        aux <- matrix(1, nrow = nrow(e.des.0.na.list[[i]]), 1)
        auxx <- matrix(1, nrow = nrow(e.des.1.list[[i]]), 1)
      }

      if (i == 1) {
        e.des.0.na <- aux
        e.des.1 <- auxx
      } else {
        e.des.0.na <- Matrix::bdiag(e.des.0.na, aux)
        e.des.1 <- Matrix::bdiag(e.des.1, auxx)
      }
    }

    rownames(e.des.0.na) <- names0
    rownames(e.des.1) <- names1
    colnames(e.des.0.na) <- colnames(e.des.1) <- paste0(names(w.constr.list), ".")

    if (verbose && (e.order > 0 || e.lags > 0)) {
      warning(paste0("One of e.order > 0 and e.lags > 0 was specified, however the current number of observations (",
                    T.e, ") used to estimate conditional moments of the out-of-sample error",
                    " is not larger than the number of parameters used in estimation (", params.e, ") plus 10.",
                    " To avoid over-fitting issues e.order and e.lags were set to 0."),
                    immediate. = TRUE, call. = FALSE)
    }
    e.order <- 0
    e.lags <- 0
  }
  
  params.e <- ncol(e.des.0.na)

  e.lb.gau <- e.ub.gau <- e.lb.ls <- e.ub.ls <- e.lb.qreg <- e.ub.qreg <- e.mean <- e.var <- c()

  e.des.0.na.list <- mat2list(e.des.0.na)
  e.des.1.list <- mat2list(e.des.1)
  e1.rnames <- rownames(e.des.1)
  e0.rnames <- rownames(e.des.0.na)
  
  for (i in seq_len(I)) {
    e.des.0.na.list[[i]] <- detectConstant(e.des.0.na.list[[i]])
    e.des.1.list[[i]] <- detectConstant(e.des.1.list[[i]], scale.x)

    if (i == 1) {
      e.des.0.na <- e.des.0.na.list[[i]]
      e.des.1 <- e.des.1.list[[i]]
    } else {
      e.des.0.na <- Matrix::bdiag(e.des.0.na, e.des.0.na.list[[i]])
      e.des.1 <- Matrix::bdiag(e.des.1, e.des.1.list[[i]])
    }
  }

  e.des.0.na <- as.matrix(e.des.0.na)
  e.des.1 <- as.matrix(e.des.1)
  rownames(e.des.0.na) <- e0.rnames
  rownames(e.des.1) <- e1.rnames

  if (e.method == "gaussian" || e.method == "all") {

    pi.e   <- scpi.out(res = e.res.na, x = e.des.0.na, eval = e.des.1,
                       e.method = "gaussian", alpha = e.alpha / 2,
                       e.lb.est = e.lb.est, e.ub.est =  e.lb.est,
                       effect = sc.pred$data$specs$effect, out.feat = out.feat)

    e.lb.gau <- pi.e$lb
    e.ub.gau <- pi.e$ub
    e.mean <- pi.e$e.1
    e.var <- pi.e$e.2

    # Overwrite with user's input
    if (e.lb.est == FALSE) e.lb.gau <- e.bounds[, 1]
    if (e.ub.est == FALSE) e.ub.gau <- e.bounds[, 2]

    sc.l.1 <- sc.l.0 + e.lb.gau - epsk
    sc.r.1 <- sc.r.0 + e.ub.gau + epsk
    len.1  <- sc.r.1 - sc.l.1
  }

  if (e.method == "ls" || e.method == "all") {

    pi.e   <- scpi.out(res = e.res.na, x = e.des.0.na, eval = e.des.1,
                       e.method = "ls", alpha = e.alpha / 2,
                       e.lb.est = e.lb.est, e.ub.est =  e.lb.est,
                       effect = sc.pred$data$specs$effect, out.feat = out.feat)

    e.lb.ls <- pi.e$lb
    e.ub.ls <- pi.e$ub
    # Overwrite with user's input
    if (e.lb.est == FALSE) e.lb.ls <- e.bounds[, 1]
    if (e.ub.est == FALSE) e.ub.ls <- e.bounds[, 2]

    sc.l.2 <- sc.l.0 + e.lb.ls - epsk
    sc.r.2 <- sc.r.0 + e.ub.ls + epsk
    len.2  <- sc.r.2 - sc.l.2

  }

  if (e.method == "qreg" || e.method == "all") {

    if (e.order == 0) {

      lb <- quantile(e.res.na, e.alpha / 2)
      ub <- quantile(e.res.na, 1-e.alpha / 2)

    } else {
      pi.e   <- scpi.out(res = e.res.na, x = e.des.0.na, eval = e.des.1,
                         e.method = "qreg", alpha = e.alpha / 2,
                         e.lb.est = e.lb.est, e.ub.est =  e.ub.est, verbose = verbose,
                         effect = sc.pred$data$specs$effect, out.feat = out.feat)

      lb <- pi.e$lb
      ub <- pi.e$ub
    }

    e.lb.qreg <- lb
    e.ub.qreg <- ub

    # Overwrite with user's input
    if (e.lb.est == FALSE) e.lb.qreg <- e.bounds[, 1]
    if (e.ub.est == FALSE) e.ub.qreg <- e.bounds[, 2]

    sc.l.3 <- sc.l.0 + e.lb.qreg - epsk
    sc.r.3 <- sc.r.0 + e.ub.qreg + epsk
    len.3  <- sc.r.3 - sc.l.3
  }

  #############################################################################
  #############################################################################
  ## Simultaneous Prediction Intervals (for each unit)
  #############################################################################
  #############################################################################
  if (sc.effect == "unit-time") { # joint within unit
    joint.bounds <- simultaneousPredGet(vsig, T1, nrow(P.na), I, u.alpha, e.alpha,
                                        e.res.na, e.des.0.na, e.des.1, w.lb.est, w.ub.est,
                                        w.bounds, w.constr.inf[[1]]["name"],
                                        sc.pred$data$specs$effect, out.feat)

  } else if (sc.effect == "unit") { # joint across units
    joint.bounds <- simultaneousPredGet(vsig, nrow(P.na), nrow(P.na), I = 1, u.alpha, e.alpha,
                                        e.res.na, e.des.0.na, e.des.1, w.lb.est, w.ub.est, w.bounds,
                                        w.constr.inf[[1]]["name"], sc.pred$data$specs$effect, out.feat)

  } else if (sc.effect == "time") { # joint within aggregate unit
    joint.bounds <- simultaneousPredGet(vsig, min(unlist(T1)), nrow(P.na), 1, u.alpha, e.alpha,
                                        e.res.na, e.des.0.na, e.des.1, w.lb.est, w.ub.est,
                                        w.bounds, w.constr.inf[[1]]["name"],
                                        sc.pred$data$specs$effect, out.feat)
  }

  ML <- joint.bounds$ML
  MU <- joint.bounds$MU

  names(ML) <- rownames(P.na)
  names(MU) <- rownames(P.na)

  if (sc.pred$data$specs$effect == "time") {
    rownames(Wlb) <- paste0("aggregate.", rownames(Wlb))
    rownames(Wub) <- paste0("aggregate.", rownames(Wub))
    names(ML) <- paste0("aggregate.", names(ML))
    names(MU) <- paste0("aggregate.", names(MU))
  }

  ## Store all bounds
  bounds <- list("insample" = cbind(Wlb, Wub),
                 "subgaussian" = cbind(Wlb + e.lb.gau - epsk, Wub + e.ub.gau + epsk),
                 "ls" = cbind(Wlb + e.lb.ls - epsk, Wub + e.ub.ls + epsk),
                 "qreg" = cbind(Wlb + e.lb.qreg - epsk, Wub + e.ub.qreg + epsk),
                 "joint" = cbind(ML - epsk.j, MU + epsk.j))

  #############################################################################
  #############################################################################
  ## Return objects
  #############################################################################
  #############################################################################

  CI.0 <- cbind(sc.l.0, sc.r.0, len.0)
  colnames(CI.0) <- c("Left Bound", "Right Bound", "Length")
  if (sc.pred$data$specs$effect == "time") rownames(P.na) <- paste0("aggregate.", rownames(P.na))
  rownames(CI.0) <- rownames(P.na)

  if (e.method == "gaussian" || e.method == "all") {
    CI.1 <- cbind(sc.l.1, sc.r.1, len.1)
    colnames(CI.1) <- c("Left Bound", "Right Bound", "Length")
    rownames(CI.1) <- rownames(P.na)
  } else {
    CI.1 <- NULL
  }

  if (e.method == "ls" || e.method == "all") {
    CI.2 <- cbind(sc.l.2, sc.r.2, len.2)
    colnames(CI.2) <- c("Left Bound", "Right Bound", "Length")
    rownames(CI.2) <- rownames(P.na)
  } else {
    CI.2 <- NULL
  }

  if (e.method == "qreg" || e.method == "all") {
    CI.3 <- cbind(sc.l.3, sc.r.3, len.3)
    colnames(CI.3) <- c("Left Bound", "Right Bound", "Length")
    rownames(CI.3) <- rownames(P.na)
  } else {
    CI.3 <- NULL
  }

  u.user <- is.null(u.design) == FALSE # True if user provided the design matrix for in-sample inference
  e.user <- is.null(e.design) == FALSE # True if user provided the design matrix for out-of-sample inference

  inference.results <- list(  CI.in.sample    = CI.0,
                              CI.all.gaussian = CI.1,
                              CI.all.ls       = CI.2,
                              CI.all.qreg     = CI.3,
                              bounds          = bounds,
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
                              rho             = rho.vec,
                              Q.star          = Q.star,
                              u.alpha         = u.alpha,
                              e.alpha         = e.alpha,
                              epskappa        = epsk,
                              sims            = sims,
                              failed.sims     = failed.sims)


  result <- list( data               = sc.pred$data,
                  est.results        = sc.pred$est.results,
                  inference.results  = inference.results)

  class(result) <- 'scpi'
  if (class.type == 'scpi_data') {
    result$data$specs$class.type <- 'scpi_scpi'
  } else if (class.type == 'scpi_data_multi') {
    result$data$specs$class.type <- 'scpi_scpi_multi'
  }

  #############################################################################
  #############################################################################
  ## Plot

  if (plot == TRUE) {
    if (is.null(plot.name) == FALSE) {
      fig.name <- plot.name
    } else {
      fig.name <- "scpi_default_plot"
    }

    if (class.type == "scpi_data") {

      scplot(result = result, fig.path = getwd(),
             fig.name = fig.name, fig.format = "png", save.data = save.data)

    } else {

      scplotMulti(result)

    }
  }

  return(result)
}
