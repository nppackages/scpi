#' @title Data Preparation for \code{scest} or \code{scpi} for Point Estimation and Inference Procedures Using Synthetic Control Methods.
#'
#' @description The command prepares the data to be used by \code{\link{scest}} or \code{\link{scpi}} to implement estimation and
#' inference procedures for Synthetic Control (SC) methods.
#' It allows the user to specify the outcome variable, the features of the treated unit to be
#' matched, and covariate-adjustment feature by feature. The names of the output matrices
#' follow the terminology proposed in \href{https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, Feng, and Titiunik (2021)}.
#'
#' Companion \href{https://www.stata.com/}{Stata} and \href{https://www.python.org/}{Python} packages are described in
#' \href{https://arxiv.org/abs/2202.05984}{Cattaneo, Feng, Palomba, and Titiunik (2022)}.
#'
#' Companion commands are: \link{scdataMulti} for data preparation in the multiple treated units case with staggered adoption,
#' \link{scest} for point estimation, \link{scpi} for inference procedures, \link{scplot} and \link{scplotMulti} for plots in
#' the single and multiple treated unit(s) cases, respectively.
#'
#' Related Stata, R, and Python packages useful for inference in SC designs are described in the following website:
#'
#' \href{https://nppackages.github.io/scpi/}{ https://nppackages.github.io/scpi/}
#'
#' For an introduction to synthetic control methods, see \href{https://www.aeaweb.org/articles?id=10.1257/jel.20191450}{Abadie (2021)} and references therein.
#'
#' @param df a dataframe object.
#' @param id.var a character or numeric scalar with the name of the variable containing units' IDs. The ID variable can be numeric or character.
#' @param time.var a character with the name of the time variable. The time variable has to be numeric, integer, or Date. In
#' case \code{time.var} is Date it should be the output of \code{\link{as.Date}()} function. An integer or
#' numeric time variable is suggested when working with yearly data, whereas for all other formats a Date type
#' time variable is preferred.
#' @param outcome.var a character with the name of the outcome variable. The outcome variable has to be numeric.
#' @param period.pre a numeric vector that identifies the pre-treatment period in time.var.
#' @param period.post a numeric vector that identifies the post-treatment period in time.var.
#' @param unit.tr a character or numeric scalar that identifies the treated unit in \code{id.var}.
#' @param unit.co a character or numeric vector that identifies the donor pool in \code{id.var}.
#' @param features a character vector containing the name of the feature variables used for estimation.
#' If this option is not specified the default is \code{features = outcome.var}.
#' @param cov.adj a list specifying the names of the covariates to be used for adjustment for each feature. If \code{outcome.var} is
#' not in the variables specified in \code{features}, we force \code{cov.adj<-NULL}. See the \strong{Details} section for more.
#' @param constant a logical which controls the inclusion of a constant term across features. The default value is \code{FALSE}.
#' @param cointegrated.data a logical that indicates if there is a belief that the data is cointegrated or not.
#' The default value is \code{FALSE}.  See the \strong{Details} section for more.
#' @param anticipation a scalar that indicates the number of periods of potential anticipation effects. Default is 0.
#' @param verbose if \code{TRUE} prints additional information in the console.
#'
#' @return
#' The command returns an object of class 'scdata' containing the following
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
#' \item{\code{outcome.var}, a character with the name of the outcome variable}
#' \item{\code{features}, a character vector with the name of the features}
#' \item{\code{constant}, for internal use only}
#' \item{\code{out.in.features}, for internal use only}
#' \item{\code{effect}, for internal use only}
#' \item{\code{sparse.matrices}, for internal use only}
#' \item{\code{treated.units}, list containing the IDs of all treated units}}}
#'
#' @details
#'
#' \itemize{
#' \item{\code{cov.adj} can be used in two ways. First, if only one feature is specified through the option \code{features},
#' \code{cov.adj} has to be a list with one (even unnamed) element (eg. \code{cov.adj = list(c("constant","trend"))}).
#' Alternatively, if multiple features are specified, then the user has two possibilities:
#' \itemize{
#' \item{provide a list with one element, then the same covariates are used for
#' adjustment for each feature. For example, if there are two features specified and the user inputs
#' \code{cov.adj = list(c("constant","trend"))}, then a constant term and a linear trend are for adjustment for both features.}
#' \item{provide a list with as many elements as the number of features specified, then
#' feature-specific covariate adjustment is implemented. For example,
#' \code{cov.adj = list('f1' = c("constant","trend"), 'f2' = c("trend"))}. In this case the name of each element
#' of the list should be one (and only one) of the features specified. Note that if two (or more) features are 
#' specified and covariates adjustment has to be specified just for one of them, the user must still provide a list 
#' of the same length of the number of features, e.g., \code{cov.adj = list('f1' = c("constant","trend"), 'f2' = NULL}.}
#' }
#'
#' This option allows the user to include feature-specific constant terms
#' or time trends by simply including "constant" or "trend" in the corresponding
#' element of the list.
#'  
#' When \code{outcome.var} is not included in \code{features}, we automatically set \eqn{\mathcal{R}=\emptyset}, that is
#' we do not perform covariate adjustment. This is because, in this setting it is natural to create the out-of-sample
#' prediction matrix \eqn{\mathbf{P}} using the post-treatment outcomes of the donor units only.
#' }
#'
#' \item{\code{cointegrated.data} allows the user to model the belief that \eqn{\mathbf{A}} and \eqn{\mathbf{B}} form a
#' cointegrated system. In practice, this implies that when dealing with the pseudo-true
#' residuals \eqn{\mathbf{u}}, the first-difference of \eqn{\mathbf{B}} are used rather than the levels.}
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
#' \item{\href{https://www.aeaweb.org/articles?id=10.1257/jel.20191450}{Abadie, A. (2021)}.
#' Using synthetic controls: Feasibility, data requirements, and methodological aspects.
#' \emph{Journal of Economic Literature}, 59(2), 391-425.}
#' \item{\href{https://nppackages.github.io/references/Cattaneo-Feng-Titiunik_2021_JASA.pdf}{Cattaneo, M. D., Feng, Y., and Titiunik, R.
#' (2021)}. Prediction intervals for synthetic control methods. \emph{Journal of the American Statistical Association}, 116(536), 1865-1880.}
#' \item{\href{https://arxiv.org/abs/2202.05984}{Cattaneo, M. D., Feng, Y., Palomba F., and Titiunik, R. (2022).}
#' scpi: Uncertainty Quantification for Synthetic Control Methods, \emph{arXiv}:2202.05984.}
#' \item{\href{https://arxiv.org/abs/2210.05026}{Cattaneo, M. D., Feng, Y., Palomba F., and Titiunik, R. (2022).}
#' Uncertainty Quantification in Synthetic Controls with Staggered Treatment Adoption, \emph{arXiv}:2210.05026.}
#'}
#'
#' @seealso \code{\link{scdataMulti}}, \code{\link{scest}}, \code{\link{scpi}}, \code{\link{scplot}}, \code{\link{scplotMulti}}
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
#' @export

scdata <- function(df,
                   id.var,
                   time.var,
                   outcome.var,
                   period.pre,
                   period.post,
                   unit.tr,
                   unit.co,
                   features = NULL,
                   cov.adj = NULL,
                   cointegrated.data = FALSE,
                   anticipation = 0,
                   constant = FALSE,
                   verbose = TRUE) {

  ############################################################################
  ############################################################################
  ### Error Checking

  # Safe copy
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


  if (is.null(cov.adj) == FALSE) {
    if (is.list(cov.adj) == FALSE) {
      stop("The argument cov.adj should be a list!")
    }

    # Check all covariates other than constant and trend are in dataframe
    unique.covs <- unique(unlist(cov.adj))
    unique.covs <- unique.covs[!unique.covs %in% c("constant", "trend")]

    if (length(unique.covs > 0)) { # not only constant and trend specified
      if (!all(unique.covs %in% var.names)) {
        cov.adj.not.found <- unique.covs[!(unique.covs %in% var.names)]
        stop(paste(c("The following covariate(s) for adjustment are not in the input dataframe:",
                     cov.adj.not.found), collapse = " "))
      }
    }

    # Check that there are no duplicates within same equation
    for (j in seq_len(length(cov.adj))) {
      cc <- cov.adj[[j]]
      aux <- table(cc)

      if (any(aux > 1)) {
        stop(paste(c("Equation ", j, " contains more than once the same covariate!")))
      }
    }
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

  time.var.class <- var.class[var.names == time.var]
  if (!time.var.class %in% c("numeric", "integer", "Date")) {
    stop("Time variable (time.var) must be either numeric or character!")
  }

  if (!var.class[var.names == outcome.var] %in% c("numeric", "integer")) {
    stop("Outcome variable (outcome.var) must be numeric!")
  }

  if (is.null(features) ==  TRUE) {
    features <- outcome.var
  }

  if (!(all(is.character(features)))) {
    stop("The object features should contain only character entries!")
  }

  if (!(all(features %in% var.names))) {
    fe.not.found <- features[!(features %in% var.names)]
    stop(paste(c("The following features are not in the input dataframe:",
                 fe.not.found), collapse = " "))
  }

  if (outcome.var %in% features) {
    out.in.features <- TRUE
  } else {
    out.in.features <- FALSE
  }

  # Error checking when the user specifies the adj eq-by-eq
  if (length(cov.adj) > 1)  {

    if (is.null(names(cov.adj))) {
      stop("You should specify the name of the feature each covariate adjustment refers to!
           (eg. cov.adj = list('feature1' = c('cov1','cov2'), 'feature2' = c('constant','cov1')))")
    }

    if (!all(names(cov.adj) %in% features)) {
      stop(paste(c("When specifying covariate adjustment separately for each feature make sure
         that there is a one-to-one match between equation names and feature names.")))
    }

    if (length(cov.adj) != length(features)) {
      stop(paste(c("When specifying covariate adjustment separately for each feature make sure
         to do it for all features! You specified covariate adjustment for ", length(cov.adj), " features
                 when you currently have ", length(features), " features!")))
    }
  }

  if ((length(features) == 1) && (constant == TRUE) && ("constant" %in% unlist(cov.adj))) {
    stop("When specifying just one feature you either specify constant == TRUE or include 'constant' in
         cov.adj!")
  }

  period.pre  <- sort(period.pre, decreasing = FALSE)
  period.post <- sort(period.post, decreasing = FALSE)

  # Create ID and time variables
  if (is.numeric(data[id.var])) {
    data[id.var] <- as.character(data[id.var])
  }

  id          <- as.matrix(unique(data[id.var]))        # ID of units
  time        <- unique(data[,time.var])                # Time periods
  time        <- time[time %in% c(period.pre, period.post)]
  
  # Check that specified units are in dataframe
  if (!(unit.tr %in% id)) {
    stop("There is no treated unit with the specified ID (unit.tr) in the specified ID variable (id.var)!")
  }

  if (!all(unit.co %in% id)) {
    co.not.found <- unit.co[!(unit.co %in% id)]

    stop(paste(c("The following control unit(s) are not in the input dataframe:",
                 co.not.found), collapse = " "))
  }

  if (unit.tr %in% unit.co) {
    stop("The treated unit is also contained in the donor pool!")
  }

  if (length(unit.co) < 2) {
    stop("Please provide at least two control units!")
  }

  # Check specified time periods are in dataframe
  if (!all(period.pre %in% time)) {
    pre.not.found <- period.pre[!(period.pre %in% time)]
    stop(paste(c("The following pre-treatment period(s) are not in the input dataframe:",
                 pre.not.found), collapse = " "))
  }

  if (!all(period.post %in% time)) {
    post.not.found <- period.post[!(period.post %in% time)]
    stop(paste(c("The following post-treatment period(s) are not in the input dataframe:",
                 post.not.found), collapse = " "))
  }

  if (any(period.pre %in% period.post)) {
    stop("There is an overlap between the pre-treatment period and post-treatment period!")
  }

  # Consider eventual anticipation effect
  if (!is.numeric(anticipation)) {
    stop("The object 'anticipation' has to be an integer!")
  }

  if (anticipation > 0) {
    t0 <- length(period.pre); d <- anticipation
    period.post <- c(period.pre[(t0-d+1):t0], period.post)
    period.pre  <- period.pre[1:(t0-d)]
  }

  # Order outcome variable as first feature and handle covariates for adjustment accordingly if needed
  if (out.in.features == TRUE) {
    if (length(cov.adj) > 1) {
      cov.adj  <- c(cov.adj[outcome.var], cov.adj[names(cov.adj) != outcome.var])
    }
    features <- c(features[features == outcome.var], features[features != outcome.var])
  }

  # Rename time and ID variables
  var.names[var.names == time.var] <- "Time"
  var.names[var.names == id.var]   <- "ID"
  names(data) <- var.names

  ############################################################################
  ############################################################################
  ### Data preparation

  # Make the panel balanced
  data.bal <- as.data.frame(tidyr::complete(data, .data[["ID"]],.data[["Time"]]))

  # Identify rows corresponding to treatment unit
  rows.tr.pre <- which(data.bal[, "ID"]   %in% c(unit.tr) &
                         data.bal[, "Time"] %in% period.pre)

  # Identify rows corresponding to control units
  rows.co.pre <- which(data.bal[, "ID"]   %in% c(unit.co) &
                         data.bal[, "Time"] %in% period.pre)

  # Identify rows corresponding to treatment unit
  rows.tr.post <- which(data.bal[, "ID"]   %in% c(unit.tr) &
                          data.bal[, "Time"] %in% period.post)

  # Identify rows corresponding to control units
  rows.co.post <- which(data.bal[, "ID"]   %in% c(unit.co) &
                          data.bal[, "Time"] %in% period.post)
  A           <- as.matrix(c(as.matrix(data.bal[rows.tr.pre, features]))) # Stack features
  colnames(A) <- unit.tr

  ### Estimation Data
  # Actual Pre-treatment Series
  Y.pre           <- as.matrix(data.bal[rows.tr.pre, outcome.var])
  rownames(Y.pre) <- paste(unit.tr,as.character(data.bal[rows.tr.pre, "Time"]),sep=".")
  colnames(Y.pre) <- outcome.var

  # Create B
  sel <- data.bal[rows.co.pre, c(features, "ID", "Time")] # Select rows and columns

  # Stack all features one on top of the other

  aux <- stats::reshape(sel,
                 direction = "long",
                 varying   = list(names(sel)[1:length(features)]),
                 idvar     = c("ID", "Time"),
                 timevar   = "feature",
                 times     = features)

  # make the df wide so that countries are one next to the other
  aux <- stats::reshape(aux,
                 direction = "wide",
                 idvar     = c("Time", "feature"),
                 timevar   = "ID")

  time.vec    <- aux[,"Time"]
  feature.vec <- aux[,"feature"]
  B           <- as.matrix(aux[, !(colnames(aux) %in% c("Time", "feature"))])

  B.names     <- stringr::str_remove(colnames(B), features[1])
  B.names     <- paste(unit.tr, stringr::str_remove(B.names,"."), sep = ".")
  colnames(B) <- B.names
  B.names     <- sort(B.names)
  B           <- B[,B.names]

  ## Create matrix with outcome for the donors
  sel <- data.bal[rows.co.pre, c(outcome.var, "ID", "Time")] # Select rows and columns
  aux <- stats::reshape(sel,
                 direction = "wide",
                 idvar     = "Time",
                 timevar   = "ID")

  # if "Time" is not in numeric format the matrix becomes a matrix of characters,
  # this is why we do everything in one step
  Y.donors <- as.matrix(aux[, colnames(aux) != "Time", drop=FALSE])

  Y.names     <- stringr::str_remove(colnames(Y.donors), outcome.var)
  Y.names     <- stringr::str_remove(Y.names,".")

  colnames(Y.donors) <- paste(unit.tr, Y.names, sep = ".")
  rownames(Y.donors) <- paste(unit.tr, as.character(aux[,'Time']), sep = ".")
  Y.donors    <- Y.donors[ , B.names]  # Re-order according to B

  ## Create C
  C        <- NULL

  if (is.null(cov.adj) == FALSE) {
    all.covs <- unique(unlist(cov.adj))  # Store all covariates used in at least one equation

    ### AAA: currently C comes from the first unit only

    # Fixed covariate adjustment across features
    if (length(cov.adj) == 1) {
      covs.adj <- unlist(cov.adj)

      # Check that constant/time trend are required by the user
      if ("constant" %in% all.covs) {
        covs.adj <- covs.adj[covs.adj != "constant"]
        C       <- cbind(C, rep(1, length(rows.tr.pre)))
      }

      if ("trend" %in% all.covs) {
        covs.adj    <- covs.adj[covs.adj != "trend"]
        time.trend  <- period.pre - period.pre[1] + 1  # It takes into account scattered data
        C           <- cbind(C, time.trend)
      }

      rows.C <- which(data.bal[, "ID"]   == unit.co[1] &
                        data.bal[, "Time"] %in% period.pre)

      C <- cbind(C, as.matrix(data.bal[rows.C, covs.adj]))

      C <- kronecker(diag(length(features)), C)


      # Feature specific covariate adjustment
    } else if (length(cov.adj) > 1) {

      for (m in seq_len(length(features))) {
        C.m          <- NULL                # Create empty block
        feature.covs <- cov.adj[[m]]        # select feature-specific list

        # Check that constant/time trend are required by the user
        if ("constant" %in% feature.covs) {
          feature.covs <- feature.covs[feature.covs != "constant"]
          C.m          <- cbind(C.m, rep(1, length(rows.tr.pre)))
        }

        if ("trend" %in% feature.covs) {
          feature.covs <- feature.covs[feature.covs != "trend"]
          time.trend   <- period.pre - period.pre[1] + 1  # It takes into account scattered data
          C.m          <- cbind(C.m, time.trend)
        }

        rows.C <- which(data.bal[, "ID"]   == unit.co[1] &
                          data.bal[, "Time"] %in% period.pre)

        C.m <- cbind(C.m, as.matrix(data.bal[rows.C, feature.covs]))

        if (m == 1) {
          C <- C.m
        } else {
          C <- Matrix::bdiag(C, C.m)
        }
      }
    }

    C <- as.matrix(C)
    colnames.C <- list()

    for (j in seq_len(length(features))) {
      if (length(cov.adj) == 1) {
        cc <- cov.adj[[1]]
      } else {
        cc           <- cov.adj[[j]]
      }

      # Check if trend and/or constant are specified and order them 2nd and 1st
      if ("trend" %in% cc) {
        cc <- c(cc[cc == "trend"], cc[cc != "trend"])
      }

      if ("constant" %in% cc) {
        cc <- c(cc[cc == "constant"], cc[cc != "constant"])
      }


      if (is.null(cc) == TRUE) {
        colnames.C[[j]] <- NULL
      } else {
        colnames.C[[j]] <- paste(unit.tr, features[j], cc, sep = ".")
      }
    }

    colnames(C) <- unlist(colnames.C)

  }

  ############################################################################
  ##############################################################################
  ### Prediction Data

  ## Actual post-treatment series
  Y.post <- as.matrix(data.bal[rows.tr.post, outcome.var])
  rownames(Y.post) <- paste(unit.tr, as.character(data.bal[rows.tr.post, "Time"]), sep = ".")
  colnames(Y.post) <- outcome.var

  ## Prediction Matrix
  # Select series of donors
  aux    <- data.bal[rows.co.post, c(outcome.var, "ID", "Time")] # Select rows and columns

  # make the df wide so that countries are one next to the other
  aux <- stats::reshape(aux,
                 direction = "wide",
                 idvar     = "Time",
                 timevar   = "ID")

  P <- as.matrix(aux[, names(aux) != "Time"])
  rownames(P) <- paste(unit.tr, as.character(aux[,'Time']), sep = ".")

  P.names     <- stringr::str_remove(colnames(P), outcome.var)
  P.names     <- stringr::str_remove(P.names,".")
  colnames(P) <- paste(unit.tr, P.names, sep = ".")
  P           <- P[,B.names, drop = F]  # Re-order as the matrix B

  # If the outcome variable is within the specified features then we need to
  # augment P with the corresponding (eventual) covariates used for adjustment.
  # If instead the outcome variable is not within the specified features
  # P is composed by the outcome variable of the donors only
  
  if (out.in.features == TRUE) {
    # Check that global constant is required by the user
    if (constant == TRUE) {
      P           <- cbind(P, rep(1, length(rows.tr.post)))
      colnames(P) <- c(colnames(P[, 1:(dim(P)[2] - 1), drop = FALSE]), paste(unit.tr,"0.constant", sep = "."))
    }
    
    # Add covariates used for adjustment in outcome variable equation (if present)
    
    if (is.null(cov.adj[[1]]) == FALSE) {
      covs.adj <- unlist(cov.adj[[1]])

      # Check that constant/time trend are required by the user
      if ("constant" %in% covs.adj) {
        covs.adj <- covs.adj[covs.adj != "constant"]
        P        <- cbind(P, rep(1, length(rows.tr.post)))
        colnames(P) <- c(colnames(P[, 1:(dim(P)[2] - 1), drop = FALSE]), paste(unit.tr,"1.constant", sep = "."))
      }

      if ("trend" %in% covs.adj) {
        covs.adj    <- covs.adj[covs.adj != "trend"]
        time.trend  <- period.post - period.pre[1] + 1  # It takes into account scattered data
        P           <- cbind(P, time.trend)
        colnames(P) <- c(colnames(P[, 1:(dim(P)[2] - 1), drop = FALSE]), paste(unit.tr,"1.trend", sep = "."))
      }

      rows.P <- which(data.bal[, "ID"] == unit.co[1] & data.bal[, "Time"] %in% period.post)

      P <- cbind(P, as.matrix(data.bal[rows.P, covs.adj]))
      if (length(covs.adj > 0)) {
        colnames(P) <- c(colnames(P[,1 : (dim(P)[2]-length(covs.adj)), drop = FALSE]), paste(unit.tr,1, covs.adj, sep = "."))
      }
    }
  } 
  
  T1 <- length(period.post)


  ############################################################################
  ############################################################################
  # Proceed cleaning missing data in the pre-treatment period
  # Check if there are annoying donors with ALL missing values in the pre-treatment period
  empty.cols <- colSums(is.na(B)) == nrow(B)
  
  if (sum(empty.cols) > 0) {
    names <- strsplit(colnames(B[,empty.cols]), "\\.")
    removed.cols <- unlist(lapply(names, "[[", 2))
    warn.text <- paste(c("The following donors have no observations in the pre-treatment period, hence they have been removed!",
                         removed.cols), collapse = " ")
    if (verbose == TRUE) {
      warning(warn.text)
    }
    dropped.co <- unit.co %in% removed.cols
    unit.co.eff <- unit.co[!dropped.co]
    
    B <- B[ , !dropped.co, drop = TRUE]
    
  } else {
    unit.co.eff <- unit.co
  }
  
  X <- cbind(A, B, C)

  rownames(X) <- paste(unit.tr, feature.vec, as.character(time.vec), sep = ".")
  select      <- rowSums(is.na(X)) == 0
  X.na        <- X[select, , drop = FALSE]
  if (nrow(X.na) == 0) {
    stop("Current specification has too many missing values and no observations are left!")
  }
  
  j1   <- dim(as.matrix(A))[2]  # Columns of A
  j2   <- j1 + dim(B)[2]        # Columns of B
  j3   <- j2
  if (is.null(C) == FALSE) j3   <- j2 + dim(C)[2]        # Columns of C
  
  A.na           <- X.na[, 1:j1, drop = FALSE]
  B.na           <- X.na[, (j1+1):j2, drop = FALSE]
  if (is.null(C) == TRUE) {
    C.na <- NULL
  } else {
    C.na           <- X.na[, (j2+1):j3, drop = FALSE]
  }
  feature.na.vec <- feature.vec[select]
  time.na.vec    <- time.vec[select]
  
  if (constant == TRUE) {
    C.na <- cbind(rep(1, dim(B.na)[1]), C.na)
    rownames(C.na) <- rownames(B.na)
    if (dim(C.na)[2] == 1){
      colnames(C.na) <- paste(unit.tr, "0.constant", sep = ".")
    } else {
      colnames(C.na) <- c(paste(unit.tr, "0.constant", sep = "."), colnames(C.na[,2:(dim(C.na)[2]), drop = FALSE]))
    }
  }

  # Store effective number of observations per feature
  xx <- as.data.frame(table(feature.na.vec))
  T0.features        <- c(xx$Freq)
  names(T0.features) <- xx$feature.na.vec
  T0.features        <- T0.features[match(names(T0.features), features)]
  
  ############################################################################
  ############################################################################
  # Store objects

  # Size of donor pool
  J  <- length(unit.co.eff)

  # Total number of covariates used for adjustment
  if (is.null(C.na) == FALSE) {
    KM <- dim(C.na)[2]
  } else {
    KM <- 0
  }

  # Number of features
  M <- length(features)

  # Vector containing number of covariates used for adjustment in each equation
  if (is.null(cov.adj) == TRUE) {
    K  <- rep(0, M)
  } else if (length(cov.adj) == 1) {
    K  <- rep(length(cov.adj[[1]]), M)
  } else if (length(cov.adj) > 1) {
    K  <- unlist(lapply(cov.adj, length))
  }

  names(K) <- features

  # Augment P with zeros for other equations' cov adj
  if (sum(K[-1]) > 0 & out.in.features == TRUE) {
    zeros        <- matrix(0, nrow = nrow(P), ncol = sum(K[-1]))
    P            <- cbind(P, zeros)
  }

  if (constant == TRUE) {
    K <- K + 1
  }

  if (out.in.features == FALSE) {
    C.na <- NULL
    KM <- 0
    nK <- names(K)
    K <- rep(0, length(K))
    names(K) <- nK
    colnames(P) <- colnames(B.na)
  } else {
    colnames(P) <- c(colnames(B.na), colnames(C.na))
  }

  specs <- list(J = J,
                K = K,
                KM = KM,
                M = M,
                I = 1,
                cointegrated.data = cointegrated.data,
                period.pre = period.pre,
                period.post = period.post,
                T0.features = T0.features,
                T1.outcome = T1,
                outcome.var = outcome.var,
                features = features,
                constant = constant,
                out.in.features = out.in.features,
                treated.units = unit.tr,
                donors.units = unit.co.eff,
                effect = "unit-time",
                sparse.matrices = FALSE,
                units.est = unit.tr)

  df.sc <-     list(A = A.na,
                    B = B.na,
                    C = C.na,
                    P = P,
               P.diff = NULL, 
                Y.pre = Y.pre,
               Y.post = Y.post,
             Y.donors = Y.donors,
                specs = specs)

  class(df.sc) <- 'scdata'

  return(df.sc = df.sc)


}
