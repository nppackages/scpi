library(scpi)

###############################################################################
###############################################################################
## Auxiliary functions for testing

test.data <- function(df,
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
## Auxiliary functions for testing

test_that("no error is returned", {
  expect_no_error(test.data(constant = TRUE, sparse.matrices=FALSE, effect="unit-time",
                  cointegrated.data = TRUE, features = list(c("gdp", "trade")),
                  cov.adj = list(c("constant", "trend"))))
  expect_no_error(test.data(constant = TRUE, sparse.matrices=FALSE, effect="time",
                  cointegrated.data = TRUE, features = list(c("gdp", "trade")),
                  cov.adj = list(c("constant", "trend"))))
  expect_no_error(test.data(constant = TRUE, sparse.matrices=FALSE, effect="unit",
                  cointegrated.data = TRUE, features = list(c("gdp", "trade")),
                  cov.adj = list(c("constant", "trend"))))
  expect_no_error(test.data(constant = TRUE, sparse.matrices=FALSE, effect="unit", post=5,
                            cointegrated.data = TRUE, features = list(c("gdp", "trade")),
                            cov.adj = list(c("constant", "trend"))))
})





