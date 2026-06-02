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

test.dataMulti <- function(effect = "unit-time",
                           post.est = 2,
                           constant = TRUE,
                           cointegrated.data = TRUE) {
  data <- scpi_germany
  data$treatment <- 0
  data[(data$country == "West Germany" & data$year >= 1991), "treatment"] <- 1
  data[(data$country == "Italy" & data$year >= 1992), "treatment"] <- 1

  scdataMulti(
    data,
    id.var = "country",
    outcome.var = "gdp",
    treatment.var = "treatment",
    time.var = "year",
    constant = constant,
    cointegrated.data = cointegrated.data,
    post.est = post.est,
    units.est = c("West Germany", "Italy"),
    effect = effect
  )
}

###############################################################################
###############################################################################


test_that("an error is returned",
          {
           test_obj <- test.data()
           xx <- matrix(0,30,30)
           expect_error(scpi(xx, verbose = F))
           expect_error(scpi(test_obj, w.constr = "ols", cores = 2, verbose = F))
           expect_error(scpi(test_obj, w.constr = list(name="wrong name"), cores = 2, verbose = F))
           expect_error(scpi(test_obj, w.constr = list(lb = - Inf, p = 2, dir = "<=", Q = 1), cores = 2, verbose = F))
           expect_error(scpi(test_obj, w.constr = list(lb = - Inf, p = "L1", dir = ">=", Q = 1), cores = 2, verbose = F))
           expect_error(scpi(test_obj, w.constr = list(p = "L1", dir = "<="), cores = 2, verbose = F))
           expect_error(scpi(test_obj, w.constr = list(lb = - Inf, dir = "<=", Q = 1), cores = 2, verbose = F))
           expect_error(scpi(test_obj, w.constr = list(lb = - Inf, dir = "<="), cores = 2, verbose = F))
           expect_error(scpi(test_obj, V = xx, cores = 2, verbose = F))
           expect_error(scpi(test_obj, e.design = xx, cores = 2, verbose = F))
           expect_error(scpi(test_obj, u.design = xx, cores = 2, verbose = F))
           expect_error(scpi(test_obj, P = xx, cores = 2, verbose = F))
          })

test_that("time-effect aggregate intervals keep row names", {
  set.seed(8894)
  test_obj <- test.dataMulti(effect = "time")

  res <- scpi(
    test_obj,
    sims = 10,
    cores = 1,
    w.constr = list(name = "simplex"),
    u.missp = TRUE,
    e.method = "gaussian",
    verbose = FALSE
  )

  ci_names <- rownames(res$inference.results$CI.in.sample)
  expect_length(ci_names, nrow(test_obj$P))
  expect_true(all(grepl("^aggregate\\.", ci_names)))
  expect_equal(rownames(res$inference.results$bounds$insample), ci_names)
})
