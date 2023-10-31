library(scpi)

###############################################################################
###############################################################################
## Auxiliary functions for testing

test.data <- function(df = NULL,
                      id.var = "country", 
                      time.var = "year", 
                      outcome.var = "gdp", 
                      period.pre = (1960:1990), 
                      period.post = (1991:1997),
                      unit.tr = "West Germany", 
                      unit.co = NULL, 
                      features = NULL, 
                      cov.adj = NULL,
                      cointegrated.data = FALSE,
                      anticipation = 0,
                      constant = TRUE) {
  
  
  if (is.null(df)) df <- scpi_germany
  if (is.null(unit.co)) unit.co <- unique(df$country)[-7]    
  
  out  <-   scdata(df = df, id.var = id.var, time.var = time.var, outcome.var = outcome.var,
                   period.pre = period.pre, period.post = period.post,
                   unit.tr = unit.tr, unit.co = unit.co, cov.adj = cov.adj, features = features,
                   constant = constant, cointegrated.data = cointegrated.data)
  
  return(out) 
}

test.dataMulti <- function(df,
                           features = NULL,
                           cov.adj = NULL,
                           cointegrated.data = FALSE,
                           post.est = NULL,
                           units.est = NULL,
                           donors.est = NULL,
                           anticipation = 0,
                           effect = "unit-time",
                           constant = FALSE,
                           verbose = FALSE,
                           sparse.matrices = FALSE) {
  data <- scpi_germany
  
  # Create a second placebo treated unit
  data$treatment <- 0
  data[(data$country == "West Germany" & data$year >= 1991), "treatment"] <- 1
  data[(data$country == "Italy" & data$year >= 1992), "treatment"] <- 1
  
  df <- scdataMulti(data, id.var = "country", outcome.var = "gdp",
                    treatment.var = "treatment", time.var = "year", constant = constant,
                    cointegrated.data = cointegrated.data, features = features,
                    cov.adj = cov.adj, verbose = verbose, post.est = post.est,
                    units.est = units.est, donors.est = donors.est, anticipation = anticipation,
                    effect = effect, sparse.matrices = sparse.matrices)
  
  return(df)
}


###############################################################################
###############################################################################

test_that("scdata: an error is returned",
          {
           test_obj <- test.data()
           xx <- matrix(0,30,30)
           expect_error(scest(xx))
           expect_error(scest(test_obj, w.constr = "ols"))
           expect_error(scest(test_obj, w.constr = list(name="wrong name")))
           expect_error(scest(test_obj, w.constr = list(lb = - Inf, p = 2, dir = "<=", Q = 1)))
           expect_error(scest(test_obj, w.constr = list(lb = - Inf, p = "L1", dir = ">=", Q = 1)))
           expect_error(scest(test_obj, w.constr = list(p = "L1", dir = "<=")))
           expect_error(scest(test_obj, w.constr = list(lb = - Inf, dir = "<=", Q = 1)))
           expect_error(scest(test_obj, w.constr = list(lb = - Inf, dir = "<=")))
           expect_error(scest(test_obj, V = xx))
          })

# test no error is returned for five main constraints
test_that("scdata: no error is returned",
          {
            test_obj <- test.data()
            expect_no_error(scest(test_obj))
            expect_no_error(scest(test_obj, w.constr=list(name="lasso")))
            expect_no_error(scest(test_obj, w.constr=list(name="ridge")))
            expect_no_error(scest(test_obj, w.constr=list(name="L1-L2")))
            expect_no_error(scest(test_obj, w.constr=list(name="ols")))
            
            # add cointegration, features, and cov.adj
            test_obj <- test.data(features=c("gdp", "trade"),
                                  cov.adj=list("gdp"=c("constant"), "trade"=c("trend")),
                                  cointegrated.data=TRUE)
            expect_no_error(scest(test_obj))
            expect_no_error(scest(test_obj, w.constr=list(name="lasso")))
            expect_no_error(scest(test_obj, w.constr=list(name="ridge")))
            expect_no_error(scest(test_obj, w.constr=list(name="L1-L2")))
            expect_no_error(scest(test_obj, w.constr=list(name="ols")))
          })

test_that("scdataMulti: no error is returned",
          {
            test_obj <- test.dataMulti()
            expect_no_error(scest(test_obj))
            expect_no_error(scest(test_obj, w.constr=list(name="lasso")))
            expect_no_error(scest(test_obj, w.constr=list(name="ridge")))
            expect_no_error(scest(test_obj, w.constr=list(name="L1-L2")))
            expect_no_error(scest(test_obj, w.constr=list(name="ols")))
            
            # add cointegration, features, and cov.adj and sparse matrices
            test_obj <- test.dataMulti(features=list(c("gdp", "trade")),
                                       cov.adj=list("gdp"=c("constant"), "trade"=c("trend")),
                                       cointegrated.data=TRUE, sparse.matrices=TRUE)
            expect_no_error(scest(test_obj))
            expect_no_error(scest(test_obj, w.constr=list(name="lasso")))
            expect_no_error(scest(test_obj, w.constr=list(name="ridge")))
            expect_no_error(scest(test_obj, w.constr=list(name="L1-L2")))
            expect_no_error(scest(test_obj, w.constr=list(name="ols")))
            
            test_obj <- test.dataMulti(features=list(c("gdp", "trade")), effect = "time",
                                       cov.adj=list("gdp"=c("constant"), "trade"=c("trend")),
                                       cointegrated.data=TRUE, sparse.matrices=TRUE)
            expect_no_error(scest(test_obj))
            expect_no_error(scest(test_obj, w.constr=list(name="lasso")))
            expect_no_error(scest(test_obj, w.constr=list(name="ridge")))
            expect_no_error(scest(test_obj, w.constr=list(name="L1-L2")))
            expect_no_error(scest(test_obj, w.constr=list(name="ols")))
            
            test_obj <- test.dataMulti(features=list(c("gdp", "trade")), effect = "unit",
                                       cov.adj=list("gdp"=c("constant"), "trade"=c("trend")),
                                       cointegrated.data=TRUE, sparse.matrices=TRUE)
            expect_no_error(scest(test_obj))
            expect_no_error(scest(test_obj, w.constr=list(name="lasso")))
            expect_no_error(scest(test_obj, w.constr=list(name="ridge")))
            expect_no_error(scest(test_obj, w.constr=list(name="L1-L2")))
            expect_no_error(scest(test_obj, w.constr=list(name="ols")))
            
            })
