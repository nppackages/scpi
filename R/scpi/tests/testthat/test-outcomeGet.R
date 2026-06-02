test_that("time-effect aggregate outcome processing accepts unnamed post fits", {
  y_pre_fit <- matrix(c(1, 2, 1.5, 2.5), ncol = 1)
  rownames(y_pre_fit) <- c("A.1", "A.2", "B.1", "B.2")

  y_post_fit <- matrix(c(3, 4), ncol = 1)

  y_df <- data.frame(
    ID = rep(c("A", "B"), each = 4),
    Time = rep(1:4, 2),
    Treatment = rep(c(0, 0, 1, 1), 2),
    Y = c(1, 2, 3, 4, 1.5, 2.5, 3.5, 4.5)
  )

  processed <- scpi:::outcomeGet(
    Y.pre.fit = y_pre_fit,
    Y.post.fit = y_post_fit,
    Y.df = y_df,
    units.est = c("A", "B"),
    treated.units = c("A", "B"),
    plot.type = "time",
    anticipation = 0,
    period.post = list(A = 3:5, B = 3:5),
    sparse.matrices = FALSE
  )

  expect_equal(unique(processed$toplot$ID), "aggregate")
  expect_equal(sum(processed$toplot$Treatment == 1), 2)
})
