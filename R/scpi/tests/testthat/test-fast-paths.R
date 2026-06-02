test_that("unconstrained diagonal-weight estimation matches weighted least squares", {
  Z <- matrix(
    c(
      1, 0,
      1, 1,
      1, 2,
      1, 3,
      1, 4
    ),
    ncol = 2,
    byrow = TRUE
  )
  colnames(Z) <- c("intercept", "trend")

  A <- matrix(c(1, 2.1, 2.9, 4.2, 5), ncol = 1)
  V <- diag(c(1, 2, 1, 2, 1))

  w.constr <- list(
    lb = -Inf,
    dir = "NULL",
    p = "no norm",
    name = "ols"
  )

  expected <- stats::lm.wfit(Z, A, w = diag(V))$coefficients
  actual <- scpi:::b.est(
    A = A,
    Z = Z,
    J = ncol(Z),
    KM = 0,
    w.constr = w.constr,
    V = V
  )

  expect_equal(actual, expected, tolerance = 1e-12, ignore_attr = TRUE)
})
