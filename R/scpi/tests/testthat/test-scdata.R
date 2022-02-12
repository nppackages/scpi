library(scpi)


# - Y in feat: T
# - Covs     : 0
# - Constant : T
# - Features : 1
# - mvs      : F

test_that("1. Simple SC",
          {
            test_obj <- test.data(constant=T)
            
            T0 <- 31; T1 <- 7; J <- 16; K <- 1; names(K) <- "gdp"; M <- 1; KM <- 1
            
            # Check matrices dimensions
            expect_equal(dim(test_obj$A), c(T0, 1))
            expect_equal(dim(test_obj$B), c(T0, J))
            expect_equal(dim(test_obj$C), c(T0, KM))
            expect_equal(dim(test_obj$P), c(T1, J + 1))
            
            
            # Check presence of constant term in C and in P (last column)
            expect_equal(sum(test_obj$C), T0)
            expect_equal(sum(test_obj$P[,ncol(test_obj$P)]), T1)
            
            # Check other data matrices
            expect_equal(dim(test_obj$Y.pre), c(T0, 1))
            expect_equal(dim(test_obj$Y.post), c(T1, 1)) 
            expect_equal(dim(test_obj$Y.donors), c(T0, J))
            
            
            # Check specifics
            expect_equal(test_obj$specs$J, J)
            expect_equal(test_obj$specs$K, K)
            expect_equal(test_obj$specs$M, M)            
          })


# - Y in feat: T
# - Covs     : T
# - Constant : T
# - Features : 1
# - mvs      : F

test_that("2. Simple SC - Covariate Adjustment",
          {
            test_obj <- test.data(constant=F, cov.adj = list(c("constant", "trend")))
            
            T0 <- 31; T1 <- 7; J <- 16; K <- 2; names(K) <- "gdp"; M <- 1; KM <- 2
            
            # Check matrices dimensions
            expect_equal(dim(test_obj$A), c(T0, 1))
            expect_equal(dim(test_obj$B), c(T0, J))
            expect_equal(dim(test_obj$C), c(T0, KM))
            expect_equal(dim(test_obj$P), c(T1, J + KM))
            
            
            # Check presence of constant term in C and in P (last column)
            expect_equal(sum(test_obj$C[,1]), T0)   # Constant term
            expect_equal(sum(test_obj$P[,(J+1)]), T1)
            
            # Check other data matrices
            expect_equal(dim(test_obj$Y.pre), c(T0, 1))
            expect_equal(dim(test_obj$Y.post), c(T1, 1)) 
            expect_equal(dim(test_obj$Y.donors), c(T0, J))
            
            
            # Check specifics
            expect_equal(test_obj$specs$J, J)
            expect_equal(test_obj$specs$K, K)
            expect_equal(test_obj$specs$M, M)            
          })


# - Y in feat: T
# - Covs     : T
# - Constant : T
# - Features : 2
# - mvs      : F

test_that("3. Multi-Feature SC - Covariate Adjustment",
          {
            test_obj <- test.data(constant = T, features = c("gdp","trade"), cov.adj = list(c("constant", "trend")))
            
            T0 <- 31; T1 <- 7; J <- 16; K <- c(3,3); names(K) <- c("gdp","trade"); M <- 2; KM <- 5;  
            
            # Check matrices dimensions
            expect_equal(dim(test_obj$A), c(T0*M, 1))
            expect_equal(dim(test_obj$B), c(T0*M, J))
            expect_equal(dim(test_obj$C), c(T0*M, KM))
            expect_equal(dim(test_obj$P), c(T1, J + KM))
            
            
            # Check presence of constant term in C and in P (last column)
            expect_equal(sum(test_obj$C[,1]), T0*M) # global constant
            expect_equal(sum(test_obj$C[,2]), T0) # constant for first feature
            expect_equal(sum(test_obj$C[,4]), T0) # constant for second feature
            expect_equal(sum(test_obj$P[,(J+1)]), T1) # global constant
            expect_equal(sum(test_obj$P[,(J+2)]), T1) # constant for first feature
            expect_equal(sum(test_obj$P[,(J+4):(J+5)]), 0) # stuff relative to second feature should be 0
            
            # Check other data matrices
            expect_equal(dim(test_obj$Y.pre), c(T0, 1))
            expect_equal(dim(test_obj$Y.post), c(T1, 1)) 
            expect_equal(dim(test_obj$Y.donors), c(T0, J))
            
            
            # Check specifics
            expect_equal(test_obj$specs$J, J)
            expect_equal(test_obj$specs$K, K)
            expect_equal(test_obj$specs$M, M)            
          })




test_that("an error is returned",
          {
           xx <- matrix(0,30,30)
           expect_error(test.data(data = xx)) 
           expect_error(test.data(id.var = "wrong variable"))
           expect_error(test.data(time.var = "wrong variable")) 
           expect_error(test.data(outcome.var = "wrong variable")) 
           expect_error(test.data(id.var = 1)) 
           expect_error(test.data(time.var = 1)) 
           expect_error(test.data(outcome.var = 1)) 
           expect_error(test.data(period.pre = "period")) 
           expect_error(test.data(period.post = "period")) 
           expect_error(test.data(period.pre = (1950:1959))) 
           expect_error(test.data(period.post = (1950:1959))) 
           expect_error(test.data(period.pre = (1959:1990))) 
           expect_error(test.data(period.post = (1991:2004))) 
           expect_error(test.data(unit.tr = 7))
           expect_error(test.data(unit.tr = "wrong country")) 
           expect_error(test.data(id.var = "index", unit.tr = "West Germany")) 
           expect_error(test.data(id.var = "index", unit.tr = 1000)) 
           expect_error(test.data(unit.co = 2)) 
           expect_error(test.data(features = "wrong feature")) 
           expect_error(test.data(features = 1)) 
           expect_error(test.data(cov.adj = "wrong")) 
           expect_error(test.data(cov.adj = c("wrong","wrong"))) 
           expect_error(test.data(cov.adj = list("constant","constant"))) 
           expect_error(test.data(features = c("gdp","trade"), cov.adj = list(gdp = 'constant', trade = 'wrong'))) 
           expect_error(test.data(features = c("gdp","trade"), cov.adj = list(gdp = 'constant', infrate = 'trend'))) 
           
          })
