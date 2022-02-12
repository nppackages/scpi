library(scpi)

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