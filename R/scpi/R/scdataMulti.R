#' @title Data preparation to use before calling \code{scest} or \code{scpi} for point estimation and inference procedures using Synthetic Control with Multiple Treated Units and Staggered Adoption.
#'
#' @description The command prepares the data to be used by \code{\link{scest}} or \code{\link{scpi}} to implement estimation and inference procedures for Synthetic Control (SC) methods
#' in the general case of multiple treated units and staggered adoption. It is a generalization of \code{\link{scdata}}, since this latter prepares
#' the data in the particular case of a single treated unit.
#' 
#' The names of the output matrices follow the terminology proposed in 
#' \href{https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, Feng, Palomba and Titiunik (2022)}  (UPDATE LINK).
#'
#' Companion \href{https://www.stata.com/}{Stata} and \href{https://www.python.org/}{Python} packages are described in \href{https://arxiv.org/abs/2202.05984}{Cattaneo, Feng, Palomba, and Titiunik (2022)}.
#'
#' Companion commands are: \link{scdata} for data preparation in the single treated unit case, \link{scest} for point estimation, \link{scpi} for inference procedures, 
#' \link{scplot} and \link{scplotMulti} for plots in the single and multiple treated unit(s) cases, respectively.
#' 
#' Related Stata, R, and Python packages useful for inference in SC designs are described in the following website:
#' 
#' \href{ https://nppackages.github.io/scpi/}{ https://nppackages.github.io/scpi/}
#' 
#' For an introduction to synthetic control methods, see \href{https://www.aeaweb.org/articles?id=10.1257/jel.20191450}{Abadie (2021)} and references therein.
#' 
#' @param df a dataframe object.
#' @param id.var a character with the name of the variable containing units' IDs. The ID variable can be numeric or character.
#' @param time.var a character with the name of the time variable. The time variable has to be numeric, integer, or Date. In 
#' case \code{time.var} is Date it should be the output of \code{\link{as.Date}()} function. An integer or 
#' numeric time variable is suggested when working with yearly data, whereas for all other formats a Date type
#' time variable is preferred.
#' @param outcome.var a character with the name of the outcome variable. The outcome variable has to be numeric. 
#' @param treatment.var a character with the name of the variable containing the treatment assignment of each unit. The referenced
#' variable has to take value 1 if the unit is treated in that period and value 0 otherwise. Please notice that, as common in the SC
#' literature, we presume that once a unit is treated it remains treated forever. If treatment.var does not comply with this requirement
#' the command would not work as expected!
#' @param features a list containing the names of the feature variables used for estimation.
#' If features is specified, then outcome.var must be included in it. If this option is not specified the
#' default is \code{features = outcome.var}. 
#' @param cov.adj a list specifying the names of the covariates to be used for adjustment for each feature.
#' @param post.est a scalar specifying the number of post-treatment periods or a list specifying the periods 
#' for which treatment effects have to be estimated for each treated unit.
#' @param units.est a list specifying the treated units for which treatment effects have to be estimated.
#' @param constant a logical which controls the inclusion of a constant term across features. The default value is \code{FALSE}. 
#' @param cointegrated.data a logical that indicates if there is a belief that the data is cointegrated or not. The default value is \code{FALSE}. 
#' @param effect a string indicating the type of treatment effect to be estimated. Options are: 'unit-time', which estimates treatment effects for each 
#' treated unit- post treatment period combination; 'unit', which estimates the treatment effect for each unit by averaging post-treatment features over time;
#' 'time', which estimates the average treatment effect on the treated at various horizons.
#' @param anticipation a scalar that indicates the number of periods of potential anticipation effects. Default is 0. 
#' @param verbose if \code{TRUE} prints additional information in the console.
#' @param sparse.matrices if \code{TRUE} all block diagonal matrices (\eqn{\mathbf{B}}, \eqn{\mathbf{C}}, and \eqn{\mathbf{P}}) 
#' are sparse matrices. This is suggested if the dimension of the dataset is large as it will likely reduce the execution time.
#' The sparse matrices will be objects of class 'dgCMatrix' or 'lgCMatrix', thus to visualize them they need to be transformed
#' in matrices, e.g. \code{View(as.matrix(B))}.
#'
#' @return
#' The command returns an object of class 'scpi_data' containing the following
#' \item{A}{a matrix containing pre-treatment features of the treated units.}
#' \item{B}{a matrix containing pre-treatment features of the control units.}
#' \item{C}{a matrix containing covariates for adjustment.}
#' \item{P}{a matrix whose rows are the vectors used to predict the out-of-sample series for the synthetic units.}
#' \item{P.diff}{for internal use only.}
#' \item{Y.df}{a dataframe containing the outcome variable for all units.}
#' \item{Y.donors}{a matrix containing the pre-treatment outcome of the control units.}
#' \item{specs}{a list containing some specifics of the data:
#' \itemize{
#' \item{\code{J}, a list containing the number of donors for each treated unit}
#' \item{\code{K}, a list containing the number of covariates used for adjustment for each feature for each treated unit}
#' \item{\code{KM}, a list containing the total number of covariates used for adjustment for each treated unit}
#' \item{\code{M}, a list containing number of features used for each treated unit}
#' \item{\code{I}, number of treated units}
#' \item{\code{KMI}, overall number of covariates used for adjustment}
#' \item{\code{period.pre}, a list containing a numeric vector with the pre-treatment period for each treated unit}
#' \item{\code{period.post}, a list containing a numeric vector with the post-treatment period for each treated unit}
#' \item{\code{T0.features}, a list containing a numeric vector with the number of periods used in estimation for each feature for each treated unit}
#' \item{\code{T1.outcome}, a list containing the number of post-treatment periods for each treated unit}
#' \item{\code{features.list}, a list containing the name of the features for each treated unit}
#' \item{\code{outcome.var}, a character containing the name of the outcome variable}
#' \item{\code{constant}, for internal use only}
#' \item{\code{effect}, for internal use only}
#' \item{\code{anticipation}, number of periods of potential anticipation effects}
#' \item{\code{out.in.features}, for internal use only}
#' \item{\code{treated.units}, list containing the IDs of all treated units}
#' \item{\code{donors.list}, list containing the IDs of the donors of each treated unit}}}
#'
#' @details
#'
#' \itemize{
#' \item{\strong{Covariate-adjustment.} See the \strong{Details} section in \code{\link{scdata}} for further information on how
#' to specify covariate-adjustment feature-by-feature.}
#' 
#' \item{\strong{Cointegration.} \code{cointegrated.data} allows the user to model the belief that \eqn{\mathbf{A}} and \eqn{\mathbf{B}} form a
#' cointegrated system. In practice, this implies that when dealing with the pseudo-true
#' residuals \eqn{\mathbf{u}}, the first-difference of \eqn{\mathbf{B}} are used rather than the levels.}
#' 
#' \item{\strong{Effect.} \code{effect} allows the user to select between two causal quantities. The default 
#' option, \code{effect = "unit-time"}, prepares the data for estimation of 
#' \deqn{\tau_{ik},\quad k\geq, i=1,\ldots,N_1,} 
#' whereas the option \code{effect = "unit"} prepares the data for estimation of 
#' \deqn{\tau_{\cdot k}=\frac{1}{N_1} \sum_{i=1}^{N_1} \tau_{i k}}
#' which is the average effect on the treated unit across multiple post-treatment periods.}
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
#'
#' @references
#'\itemize{
#' \item{\href{https://www.aeaweb.org/articles?id=10.1257/jel.20191450}{Abadie, A. (2021)}. Using synthetic controls: Feasibility, data requirements, and methodological aspects.
#' \emph{Journal of Economic Literature}, 59(2), 391-425.}
#' \item{\href{https://cattaneo.princeton.edu/papers/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, M. D., Feng, Y., and Titiunik, R. 
#' (2021)}. Prediction intervals for synthetic control methods. \emph{Journal of the American Statistical Association}, 116(536), 1865-1880.}
#' \item{\href{https://arxiv.org/abs/2202.05984}{Cattaneo, M. D., Feng, Y., Palomba F., and Titiunik, R. (2022).}
#' scpi: Uncertainty Quantification for Synthetic Control Methods, \emph{arXiv}:2202.05984.}
#'}
#'
#' @seealso \code{\link{scdata}}, \code{\link{scest}}, \code{\link{scpi}}, \code{\link{scplot}}, \code{\link{scplotMulti}}
#'
#' @examples
#' 
#' datager <- scpi_germany
#'
#' datager$tr_id <- 0
#' datager$tr_id[(datager$country == "West Germany" & datager$year > 1990)] <- 1
#' datager$tr_id[(datager$country == "Italy" & datager$year > 1992)] <- 0
#' 
#' outcome.var <- "gdp"
#' id.var <- "country"
#' treatment.var <- "tr_id"
#' time.var <- "year"
#' df.unit <- scdataMulti(datager, id.var = id.var, outcome.var = outcome.var, 
#'                        treatment.var = treatment.var,
#'                        time.var = time.var, features = list(c("gdp", "trade")),
#'                		    cointegrated.data = TRUE, constant = TRUE)
#'
#' @export

scdataMulti <- function(df, 
                        id.var, 
                        time.var, 
                        outcome.var,
                        treatment.var,
                        features = NULL, 
                        cov.adj = NULL,
                        cointegrated.data = FALSE,
                        post.est = NULL,
                        units.est = NULL,
                        anticipation = 0,
                        effect = "unit-time",
                        constant = FALSE, 
                        verbose = TRUE,
                        sparse.matrices = FALSE) {

  ############################################################################
  ############################################################################
  ### Error Checking
  ID <- Treatment <- Time <- NULL
  
  data <- df

  # Store variable names and variable class
  var.names   <- names(data)              # Var names in dataframe
  var.class   <- sapply(data, class)      # Var types in dataframe

  # Check main input is a dataframe
  if (is.data.frame(data) == FALSE) {
    stop("Data input should be a dataframe object!")
  }

  # Convert from tibble to proper dataframe object
  if (tibble::is_tibble(data) == TRUE) {
    data <- as.data.frame(data)
  }

  # Check inputs are string
  if (is.character(id.var) == FALSE) {
    stop("You should specify the name of id.var as a character! (eg. id.var = 'ID')")
  }

  if (is.character(outcome.var) == FALSE) {
    stop("You should specify the name of outcome.var as a character! (eg. outcome.var = 'outcome')")
  }

  if (is.character(time.var) == FALSE) {
    stop("You should specify the name of time.var as a character! (eg. time.var = 'time')")
  }

  if (is.character(treatment.var) == FALSE) {
    stop("You should specify the name of treatment.var as a character! (eg. treatment.var = 'treatment')")
  }

  if (is.null(cov.adj) == FALSE){
    if (is.list(cov.adj) == FALSE) {
      stop("The argument cov.adj should be a list!")
    }
  }
  
  if (is.null(features) == FALSE) {
    if (is.list(features) == FALSE) {
      stop("The argument features should be a list!")
    }
    if (!all(names(features) %in% units.est)) {
      stop("The object features should be a list whose elements have the same name of the selected treated units!")
    }
  } else {
    features <- outcome.var
  }

  # Check inputs are in dataframe
  if (!(id.var %in% var.names)) {
    stop("ID variable (id.var) not found in the input dataframe!")
  }

  if (!(time.var %in% var.names)) {
    stop("Time variable (time.var) not found in the input dataframe!")
  }

  if (!(outcome.var %in% var.names)) {
    stop("Outcome variable (outcome.var) not found in the input dataframe!")
  }

  if (!(treatment.var %in% var.names)) {
    stop("Treatment variable (treatment.var) not found in the input dataframe!")
  }

  time.var.class = var.class[var.names == time.var]
  if (!time.var.class %in% c("numeric","integer","Date")) {
    stop("Time variable (time.var) must be either numeric or Date format!")
  }

  if (!var.class[var.names == outcome.var] %in% c("numeric","integer")) {
    stop("Outcome variable (outcome.var) must be numeric!")
  }

  if (!var.class[var.names == treatment.var] %in% c("numeric","integer")) {
    stop("Outcome variable (treatment.var) must be numeric!")
  }  

  if (is.character(effect) == FALSE) {
    stop("The option 'effect' must be a character (eg. effect = 'unit-time')!")
  }

  if (!(effect %in% c('unit-time', 'unit', 'time'))) {
    stop("The option 'effect' should be either 'unit-time', 'time', or 'unit'." )
  }

  # Rename treatment, time, and ID variables
  var.names[var.names == treatment.var] <- "Treatment"
  var.names[var.names == id.var]   <- "ID"
  var.names[var.names == time.var]   <- "Time"
  names(data) <- var.names
  if (is.numeric(data[,"Treatment"])) {
    data[,"ID"] <- as.character(data[,"ID"])
  }
  time <- unique(data[,"Time"])                
  Y.df <- data[c("ID", "Time", "Treatment", outcome.var)]

  # Identify treated units
  treated.count <- aggregate(data[,"Treatment"], by=list(unit=data[,"ID"]), FUN=sum)
  treated.units <- treated.count$unit[treated.count$x > 0]
  treated.post <- treated.units

  if (is.null(units.est) == FALSE) {
    if (!all(units.est %in% treated.units)) {
      stop("The object units.est must contain the identifiers (id.var) of the treated units for which treatment effects have to be estimated!")
    }

    sel.tr <- treated.units %in% units.est
    treated.units <- treated.units[sel.tr]
  }  else {
    units.est <- treated.units
  }

  # Control covariates for adjustment and matching features 
  if (is.null(cov.adj) == FALSE) {
    cov.adj.list <- list()
    for (i in treated.units) {
      cov.adj.list[[i]] <- cov.adj
    }
  } else {
    cov.adj.list <- rep(list(NULL), length(treated.units))
    names(cov.adj.list) <- treated.units
  }
  if (length(features) == 1) {
    features.list <- rep(features, length(treated.units))
    names(features.list) <- treated.units
  } else {
    features.list <- features
  }

  # Control other boolean options and anticipation effect
  constant.list <- rep(constant, length(treated.units))
  names(constant.list) <- treated.units
    
  cointegrated.data.list <- rep(cointegrated.data, length(treated.units))
  names(cointegrated.data.list) <- treated.units
    
  anticipation.list <- rep(anticipation, length(treated.units))
  names(anticipation.list) <- treated.units
  
  if (is.null(post.est) == FALSE) {
    if (is.list(post.est) == FALSE & is.numeric(post.est) == FALSE) {
      stop("The object post.est has to be a list or an integer!")
    }
    
    if (is.list(post.est) == TRUE) {
      if (length(post.est) != length(treated.units)) {
        stop(paste0("The object post.est must have the same number of elements as the number of treated units: ",
                    length(treated.units)))
      }
      if (!all(names(post.est) %in% treated.units)) {
        stop("The elements of post.est must be named with the identifier of the treated units in id.var!")
      }
    } 
  }

  
  ############################################################################
  ############################################################################
  ### Data preparation
  
  # For each treated unit identify the pre- and post-treatment periods and the donor pool
  treated.periods <- subset(data, Treatment == 1, select=c(Time, ID)) # post treatment period for each treated unit
  
  J.list <- list()
  K.list <- list()
  KM.list <- list()
  M.list <- list()
  period.pre.list <- list()
  period.post.list <- list()
  T0.features.list <- list()
  T1.list <- list()
  out.in.features.list <- list()
  donors.list <- list()
  
  rownames.A <- c()
  colnames.B <- c()
  colnames.C <- c()
  rownames.P <- c()
  colnames.P <- c()
  rownames.Y.donors <- c()
  colnames.Y.donors <- c()
  
  tr.count <- 1
  for (treated.unit in treated.units) {
    treated.unit <- as.character(treated.unit)  # handle numeric identifier for treated units
    treated.unit.T0 <- min(treated.periods$Time[treated.periods$ID == treated.unit])   # get first treatment period of treated unit
    #donors <- subset(data, ID != treated.unit & Time < treated.unit.T0, select=c(ID, Treatment))   # get all other units before treatment of treated unit
    donors <- subset(data, !(ID %in% treated.post) & Time < treated.unit.T0, select=c(ID, Treatment))   # get all other units before treatment of treated unit
    if (nrow(donors) == 0) stop(paste0("The current specification for ", treated.unit," does not have observations!"))
    donors.count <- aggregate(donors[,"Treatment"], by=list(unit=donors$ID), FUN=sum) # number of periods units have been treated before treatment of treated unit
    donors.units <- donors.count$unit[donors.count$x == 0] # donors are those unit that have never been treated before treatment of treated unit

    if (is.null(post.est) == FALSE) {
      if (is.list(post.est)) {
        T1.last <- post.est[[treated.unit]][length(post.est[[treated.unit]])]
      } else {
        T1.last <- treated.unit.T0 + post.est
      }
      
      treated.donors <- subset(data, ID %in% treated.post & Time < T1.last, select=c(ID, Treatment))   # get all other units before treatment of treated unit
      tr.donors.count <- aggregate(treated.donors[,"Treatment"], by=list(unit=treated.donors$ID), FUN=sum) # number of periods units have been treated before treatment of treated unit
      tr.donors.units <- tr.donors.count$unit[tr.donors.count$x == 0] # donors are those unit that have never been treated before treatment of treated unit
      
      donors.units <- sort(c(donors.units, tr.donors.units)) 
    }    
    
    df.aux <- subset(data, ID %in% c(donors.units, treated.unit)) # subset dataset
    time.vec <- unique(data["Time"])
    
    period.pre <- time.vec[time.vec < treated.unit.T0]
    period.post <- time.vec[time.vec >= treated.unit.T0]
    
    if (is.null(post.est) == FALSE) {
      if (is.list(post.est)) {
        sel.post <- period.post %in% post.est[[treated.unit]] 
      } else {
        sel.post <- period.post < T1.last
      }
      period.post <- period.post[sel.post]
    }

    scdata.out <- tryCatch(
                      {
                         scdata(df.aux, id.var = "ID",
                         time.var = "Time",
                         outcome.var = outcome.var,
                         period.pre = period.pre,
                         period.post = period.post,
                         unit.tr = treated.unit,
                         unit.co = donors.units,
                         cov.adj = cov.adj.list[[treated.unit]],
                         features = features.list[[treated.unit]],
                         constant = constant.list[[treated.unit]],
                         cointegrated.data = cointegrated.data.list[[treated.unit]],
                         anticipation = anticipation.list[[treated.unit]], verbose = verbose)
                      },

                      warning = function(war) {
                        message(paste("There is a warning related to your specification for the treated unit:", treated.unit))
                        message("Here's the original warning message:")
                        if (verbose) warning(war, call. = FALSE, immediate. = TRUE)

                        aux <- scdata(df.aux, id.var = "ID", 
                                      time.var = "Time", 
                                      outcome.var = outcome.var, 
                                      period.pre = period.pre, 
                                      period.post = period.post, 
                                      unit.tr = treated.unit, 
                                      unit.co = donors.units, 
                                      cov.adj = cov.adj.list[[treated.unit]], 
                                      features = features.list[[treated.unit]],
                                      constant = constant.list[[treated.unit]], 
                                      cointegrated.data = cointegrated.data.list[[treated.unit]],
                                      anticipation = anticipation.list[[treated.unit]], verbose = FALSE)
                        return(aux)
                      },

                      error = function(err) {
                        message(paste("There is a problem with your specification for the treated unit:", treated.unit))
                        message("Here's the original error message:")
                        stop(err)
                      },

                      finally = {}
                      )

    # Store matrices with data
    A.tr <- scdata.out$A
    B.tr <- scdata.out$B
    C.tr <- scdata.out$C

    if (is.null(C.tr)) { # to avoid bdiag to crash when C is null
      C.tr <- as.matrix(data.frame()[1:nrow(B.tr), ])
    }
    Y.donors.tr <- scdata.out$Y.donors
    P.tr <- scdata.out$P 
    if (!(effect == "unit" && cointegrated.data.list[[treated.unit]])) {
      P.diff <- NULL
    }

    if (effect == "unit") { # average within unit
      if (cointegrated.data.list[[treated.unit]] ==  TRUE) { # differentiate the data if cointegration

        if (scdata.out$specs$out.in.features == FALSE) {
          P.first <- P.tr[1, , drop = FALSE] - Y.donors.tr[scdata.out$specs$T0.features[1], , drop = FALSE]
          P.diff  <- rbind(P.first, apply(P.tr, 2, diff))
          P.diff <- t(as.matrix(colMeans(P.diff)))
        } else {
          feature.id <- unlist(purrr::map(stringr::str_split(rownames(B.tr), "\\."), 2))
          T0.sel <- scdata.out$specs$T0.features[outcome.var]
          JJ <- scdata.out$specs$J
          ## Take first differences of P
          # Remove last observation of first feature from first period of P
          P.first <- c((P.tr[1, (1:JJ), drop = FALSE] - B.tr[feature.id == outcome.var, , drop = FALSE][T0.sel, ]),
                      P.tr[1 ,-(1:JJ), drop = FALSE])
          # Take differences of other periods
          if (nrow(P.tr) > 2) {
            Pdiff <- apply(P.tr[, (1:JJ), drop = FALSE], 2, diff)
            P.diff <- rbind(P.first, cbind(Pdiff, P.tr[-1, -(1:JJ), drop = FALSE]))
          } else if (nrow(P.tr) == 2) {
            Pdiff <- t(as.matrix(apply(P.tr[, (1:JJ), drop = FALSE], 2, diff)))
            P.diff <- rbind(P.first, cbind(Pdiff, P.tr[-1, -(1:JJ), drop = FALSE]))
          } else {
            P.diff <- t(as.matrix(P.first))
          }
          P.diff <- t(as.matrix(colMeans(P.diff)))
        }
      }

      P.tr <- t(as.matrix(colMeans(P.tr)))
      rownames(P.tr) <- paste(treated.unit,
                              scdata.out$specs$period.post[ceiling(scdata.out$specs$T1.outcome / 2)], sep = ".")

    }

    # Handle naming of rows
    rownames.A <- c(rownames.A, rownames(A.tr))
    rownames.P <- c(rownames.P, rownames(P.tr))
    colnames.B <- c(colnames.B, colnames(B.tr))

    if (!is.null(cov.adj.list[[treated.unit]]) || constant.list[[treated.unit]]) {
      colnames.C <- c(colnames.C, colnames(C.tr))
    }
    colnames.P <- c(colnames.P, c(colnames(B.tr), colnames(C.tr)))
    colnames.Y.donors <- c(colnames.Y.donors, colnames(Y.donors.tr))
    rownames.Y.donors <- c(rownames.Y.donors, rownames(Y.donors.tr))
    
    if (effect == "time") {
      P.tr <- cbind(c(1:nrow(P.tr)), P.tr)
      colnames(P.tr) <- c("aux_id", colnames(P.tr)[-1])
    }
    
    if (tr.count == 1) {
      A.stacked <- A.tr
      B.stacked <- B.tr
      C.stacked <- C.tr
      if (effect == "time") {
        P.stacked <- as.data.frame(P.tr)
      } else {
        P.stacked <- P.tr
      }
      
      Pd.stacked <- P.diff
      Y.donors.stacked <- Y.donors.tr

    } else {
      A.stacked <- rbind(A.stacked, A.tr)
      B.stacked <- Matrix::bdiag(B.stacked, B.tr)
      C.stacked <- Matrix::bdiag(C.stacked, C.tr)
      if (effect == "time") {
        P.stacked <- merge(P.stacked, as.data.frame(P.tr), by = "aux_id", all = TRUE)
      } else {
        P.stacked <- Matrix::bdiag(P.stacked, P.tr)
      }
      if (is.null(Pd.stacked) == FALSE) Pd.stacked <- Matrix::bdiag(Pd.stacked, P.diff)
      Y.donors.stacked <- Matrix::bdiag(Y.donors.stacked, Y.donors.tr)
    }

    J.list[[treated.unit]] <- scdata.out$specs$J
    K.list[[treated.unit]] <- scdata.out$specs$K
    KM.list[[treated.unit]] <- scdata.out$specs$KM
    M.list[[treated.unit]] <- scdata.out$specs$M
    period.pre.list[[treated.unit]] <- scdata.out$specs$period.pre
    period.post.list[[treated.unit]] <- scdata.out$specs$period.post
    T0.features.list[[treated.unit]] <- scdata.out$specs$T0.features
    T1.list[[treated.unit]] <- scdata.out$specs$T1.outcome
    out.in.features.list[[treated.unit]] <- scdata.out$specs$out.in.features
    donors.list[[treated.unit]] <- scdata.out$specs$donors.units
    tr.count <- tr.count + 1
  }
  
  # Handle names of stacked matrices
  rownames(A.stacked) <- rownames.A
  rownames(B.stacked) <- rownames.A
  rownames(C.stacked) <- rownames.A
  rownames(Y.donors.stacked) <- rownames.Y.donors
  if (effect == "time") {
    P.stacked$aux_id <- NULL
  } else {
    rownames(P.stacked) <- rownames.P
  }
  
  colnames(A.stacked) <- "A"
  colnames(B.stacked) <- colnames.B
  colnames(C.stacked) <- colnames.C
  colnames(P.stacked) <- colnames.P
  colnames(Y.donors.stacked) <- colnames.Y.donors

  
  # Rearrange P so that order of columns coincides with (B,C)
  P.stacked <- P.stacked[, c(colnames.B, colnames.C), drop = FALSE]
  if (is.null(Pd.stacked) == FALSE) {
    colnames(Pd.stacked) <- colnames.P
    rownames(Pd.stacked) <- rownames.P
    Pd.stacked <- Pd.stacked[, c(colnames.B, colnames.C), drop = FALSE]
  }
  
  if (effect == "time") {
    P.stacked <- P.stacked/length(treated.units)
    P.stacked <- na.omit(P.stacked)
  }

  specs <- list(J = J.list,
                K = K.list,
                KM = KM.list,
                M = M.list,
                I = length(treated.units),
                KMI = ncol(C.stacked),
                cointegrated.data = cointegrated.data.list,
                period.pre = period.pre.list,
                period.post = period.post.list,
                T0.features = T0.features.list,
                T1.outcome = T1.list,
                outcome.var = outcome.var,
                features = features.list,
                constant = constant.list,
                out.in.features = out.in.features.list,
                treated.units = treated.units,
                donors.list = donors.list,
                effect = effect,
                anticipation = anticipation,
                units.est = units.est)

  if (sparse.matrices == FALSE) {
    if (is.null(Pd.stacked == FALSE)) Pd.stacked <- as.matrix(Pd.stacked)

    df.sc <-     list(A = A.stacked,
                      B = as.matrix(B.stacked),
                      C = as.matrix(C.stacked),
                      P = as.matrix(P.stacked),
                      P.diff = Pd.stacked,
                      Y.df = Y.df,
                      Y.donors = as.matrix(Y.donors.stacked),
                      specs = specs)

  } else {
    df.sc <-     list(A = A.stacked,
                      B = B.stacked,
                      C = C.stacked,
                      P = P.stacked,
                      P.diff = Pd.stacked,
                      Y.df = Y.df,
                      Y.donors = Y.donors.stacked,
                      specs = specs)
  }

  class(df.sc) <- 'scpi_data_multi'

  return(df.sc = df.sc)
}