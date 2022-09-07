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
                      constant = TRUE, 
                      mv = F) {
  
  
  if (is.null(df)) df <- scpi_germany
  if (is.null(unit.co)) unit.co <- unique(df$country)[-7]    
  
  out  <-   scdata(df = df, id.var = id.var, time.var = time.var, outcome.var = outcome.var,
                   period.pre = period.pre, period.post = period.post,
                   unit.tr = unit.tr, unit.co = unit.co, cov.adj = cov.adj, features = features,
                   constant = constant, cointegrated.data = cointegrated.data)
  
  return(out) 
}

###############################################################################
###############################################################################


test_that("an error is returned",
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