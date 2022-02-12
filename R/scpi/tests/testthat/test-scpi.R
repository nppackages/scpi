library(scpi)

test_that("an error is returned",
          {
           test_obj <- test.data()
           xx <- matrix(0,30,30)
           expect_error(scpi(xx))
           expect_error(scpi(test_obj, w.constr = "ols", cores = 2))
           expect_error(scpi(test_obj, w.constr = list(name="wrong name"), cores = 2))
           expect_error(scpi(test_obj, w.constr = list(lb = - Inf, p = 2, dir = "<=", Q = 1), cores = 2))
           expect_error(scpi(test_obj, w.constr = list(lb = - Inf, p = "L1", dir = ">=", Q = 1), cores = 2))
           expect_error(scpi(test_obj, w.constr = list(p = "L1", dir = "<="), cores = 2))
           expect_error(scpi(test_obj, w.constr = list(lb = - Inf, dir = "<=", Q = 1), cores = 2))
           expect_error(scpi(test_obj, w.constr = list(lb = - Inf, dir = "<="), cores = 2))
           expect_error(scpi(test_obj, V = xx, cores = 2))
           expect_error(scpi(test_obj, e.design = xx, cores = 2))
           expect_error(scpi(test_obj, u.design = xx, cores = 2))
           expect_error(scpi(test_obj, P = xx, cores = 2))
          })