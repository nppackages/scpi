test_that("diagonal in-sample uncertainty is stable across core counts", {
  skip_on_cran()
  skip_if(
    parallel::detectCores(logical = TRUE) < 2,
    "Parallel regression check requires at least two cores."
  )

  data <- scpi_germany
  unit.co <- setdiff(unique(data$country), "West Germany")
  df <- scdata(
    df = data,
    id.var = "country",
    time.var = "year",
    outcome.var = "gdp",
    period.pre = 1960:1990,
    period.post = 1991:1997,
    unit.tr = "West Germany",
    unit.co = unit.co,
    constant = FALSE,
    cointegrated.data = TRUE
  )

  set.seed(8894)
  serial <- scpi(
    df,
    sims = 100,
    w.constr = list(name = "simplex"),
    cores = 1,
    u.order = 1,
    u.lags = 0,
    u.sigma = "HC1",
    u.missp = TRUE,
    e.order = 1,
    e.lags = 0,
    e.method = "gaussian",
    verbose = FALSE
  )

  set.seed(8894)
  parallel <- scpi(
    df,
    sims = 100,
    w.constr = list(name = "simplex"),
    cores = 2,
    u.order = 1,
    u.lags = 0,
    u.sigma = "HC1",
    u.missp = TRUE,
    e.order = 1,
    e.lags = 0,
    e.method = "gaussian",
    verbose = FALSE
  )

  serial.inf <- serial$inference.results
  parallel.inf <- parallel$inference.results

  expect_equal(
    serial.inf$CI.in.sample,
    parallel.inf$CI.in.sample,
    tolerance = 1e-8
  )
  expect_equal(
    serial.inf$CI.all.gaussian,
    parallel.inf$CI.all.gaussian,
    tolerance = 1e-8
  )
  expect_equal(
    serial.inf$bounds$insample,
    parallel.inf$bounds$insample,
    tolerance = 1e-8
  )
  expect_equal(serial.inf$failed.sims, parallel.inf$failed.sims)
})
