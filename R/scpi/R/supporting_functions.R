###############################################################################
### Auxiliary functions for estimation

# Quadratic loss function
obj.fun.est <- function(x, Z, V, A, J, QQ, KM, p) {
  f <- x %*% t(Z) %*% V %*% Z %*% x - 2 * t(A) %*% V %*% Z %*% x
  g <- 2*t(Z) %*% V %*% Z %*% x - 2 * t(t(A) %*% V %*% Z)
  
  return(list("objective" = f,
              "gradient"  = g))
}

# Set constraints
norm.co.est <- function(x, Z, V, A, J, QQ, KM, p) {

  if (p == 1) {
    av <- rep(1,J)
    av[x[1:J] < 0] <- -1
    ja <- c(av, rep(0,KM))
    co <- sum(abs(x[1:J])) - QQ
  } else {
    ja <- c(2*x[1:J], rep(0,KM))
    co <- sum(x[1:J]^2) - QQ^2
  }
  
  return(list("constraints" = co,
              "jacobian"    = ja))
}

### Auxiliary functions for inference

# Prepare objective functions
obj.fun.min <- function(x, xt, beta, Q, G, J, KM, QQ, p.int) {
  f <- -sum(xt*(x - beta))
  g <- -xt
  
  return(list("objective" = f,
              "gradient"  = g))
}

obj.fun.max <- function(x, xt, beta, Q, G, J, KM, QQ, p.int) {
  f <- sum(xt*(x - beta))
  g <- xt
  
  return(list("objective" = f,
              "gradient"  = g))
}


# Prepare inequality constraint(s): loss functions constraint + (inequality norm constraint)

# Unique inequality constraint
single.ineq <- function(x, xt, beta, Q, G, J, KM, QQ, p.int) {
  a <- -2*G - 2*c(t(beta) %*% Q)
  d <- 2*sum(G*beta) + sum(beta*(Q %*% beta))
  
  co <- x %*% Q %*% x + sum(a*x) + d
  ja <- 2*Q %*% x + a
  
  return(list("constraints" = co,
              "jacobian"    = ja))
}


# Inequality constraints: loss function + L1/L2 norm
double.ineq <- function(x, xt, beta, Q, G, J, KM, QQ, p.int) {
  # Loss function constraint
  a <- -2*G - 2*c(t(beta) %*% Q)
  d <- 2*sum(G*beta) + sum(beta*(Q %*% beta))
  
  co1 <- x %*% Q %*% x + sum(a*x) + d
  ja1 <- 2*Q %*% x + a
  
  # Norm constraint
  if (p.int == 1) {
    co2 <- sum(abs(x[1:J])) - QQ
    av <- rep(1,J)
    av[x[1:J] < 0] <- -1
    ja2 <- c(av, rep(0,KM))
  } else {
    co2 <- sum(x[1:J]^2) - QQ^2
    ja2 <- c(2*x[1:J], rep(0,KM))
  }
  
  ja <- c(t(cbind(ja1, ja2)))    # vectorize matrix row-by-row  
  
  return(list("constraints" = c(co1,co2),
              "jacobian"    = ja))
}
    
# Eventual equality constraint on norm
norm.equal <- function(x, xt, beta, Q, G, J, KM, QQ, p.int) {
  
  if (p.int == 1) {
    co <- sum(abs(x[1:J])) - QQ
    av <- rep(1,J)
    av[x[1:J] < 0] <- -1
    ja <- c(av, rep(0,KM))
  } else {
    co <- sum(x[1:J]^2) - QQ^2
    ja <- c(2*x[1:J], rep(0,KM))
  }
  
  return(list("constraints" = co,
              "jacobian"    = ja))
}


# Auxiliary function that creates the constraints to be passed to the optimization problem
w.constr.OBJ <- function(w.constr, A, Z, V, J, KM) {
  # Default method to estimate weights as in Abadie et al. (2010)
  if (is.null(w.constr)) {
    w.constr <- list(lb   = 0,
                     p    = "L1",
                     dir  = "==",
                     Q    = 1,
                     name = "simplex")
    
  } else if (w.constr[["name"]] == "simplex") {
    
    if (!"Q" %in% names(w.constr)) {
      Q <- 1
    } else {
      Q <- w.constr[["Q"]]
    }
    
    w.constr <- list(lb   = 0,
                     p    = "L1",
                     dir  = "==",
                     Q    = Q,
                     name = 'simplex')
    
  } else if (w.constr[["name"]] == "ols") {
    w.constr <- list(lb   = -Inf,
                     dir  = "NULL",
                     p    = "no norm",
                     name = 'ols')
    
  } else if (w.constr[["name"]] == "lasso") {
    
    if (!"Q" %in% names(w.constr)) {
      w.constr[["Q"]] <- shrinkage.EST("lasso", A, Z, V, J, KM)$Q
    }
    
    w.constr <- list(lb   = -Inf,
                     p    = "L1",
                     dir  = "<=",
                     Q    = w.constr[["Q"]],
                     name = 'lasso')
    
  } else if (w.constr[["name"]] == "ridge") {
    
    if (!"Q" %in% names(w.constr)) {
      aux <- shrinkage.EST("ridge", A, Z, V, J, KM)
      w.constr[["Q"]]      <- aux$Q
      w.constr[["lambda"]] <- aux$lambda
    }
    
    w.constr <- list(lb     = -Inf,
                     p      = "L2",
                     dir    = "<=",
                     Q      = w.constr[["Q"]],
                     name   = 'ridge',
                     lambda = w.constr[["lambda"]])
    
  } else {
    # if constraint is entirely user specified just check everything is fine
    if (!(all(c('p','dir','Q','lb') %in% names(w.constr)))) {
      stop("If 'name' is not specified, w.constr should be a list whose elements 
            must be named 'p','dir','Q','lb'.")
    }
    
    if (!(w.constr[["p"]] %in% c("no norm","L1","L2"))) {
      stop("Specify either p = 'no norm' (no constraint on the norm of weights),
          p = 'L1' (L1-norm), p = 'L2' (L2-norm)")
    } else if (w.constr[["p"]] == "no norm") {
      w.constr[["dir"]] <- "NULL"
    }
    
    if (!(w.constr[["dir"]] %in% c("<=","==","NULL"))) {
      stop("Specify either dir = '<=' (inequality constraint on the norm of the weights)
            or dir = '==' (equality constraint on the norm of the weights) or dir = 'NULL'
            in case you don't want to specify a constraint on the norm of the weights.")
    }

    if (!(w.constr[["lb"]] == 0 | w.constr[["lb"]] == -Inf)) {
      stop("Specify either lb = 0 or lb = -Inf.")
    } 
    
    w.constr[["lb"]]   <- w.constr[["lb"]]
    w.constr[["name"]] <- "user provided"
  }
  
  return(w.constr)
}

shrinkage.EST <- function(method, A, Z, V, J, KM) {
  
  lambd <- NULL
  if (method == "lasso") Q <- 1
  
  if (method == "ridge") {
    wls     <- lm(A ~ Z - 1, weights = diag(V))
    sig.wls <- sigma(wls)
    lambd   <- sig.wls^2*(J+KM)/sqrt(sum(wls$coef^2))     # rule of thumb for lambda (Hoerl et al, 1975)
    Q       <- sqrt(sum(wls$coef^2))/(1+lambd)            # convert lambda into Q
  }
  
  return(list(Q = Q, lambda = lambd))
}


# Auxiliary function that solves the (un)constrained problem to estimate b
# depending on the desired method
b.est <- function(A, Z, J, KM, w.constr, V, QQ, opt.list) {

  dire <- w.constr[["dir"]]
  lb   <- w.constr[["lb"]]
  p    <- w.constr[["p"]]
  
  if (w.constr[["p"]] == "no norm") {pp <- 0}
  if (w.constr[["p"]] == "L1") {pp <- 1}
  if (w.constr[["p"]] == "L2") {pp <- 2}
  
  opt.list <- prepareOptions(opt.list, p, dire, lb, "scest")
  
  use.CVXR <- useCVXR(w.constr)
  
  if (use.CVXR == TRUE) { # handle L1 norm + inequality constraint
    
    x <- CVXR::Variable(J+KM)

    omega.inv <- base::solve(V)
    C         <- chol(omega.inv)
    C.inv     <- base::solve(C)
    A.star    <- C.inv %*% A
    Z.star    <- C.inv %*% Z
    
    objective   <- CVXR::Minimize(CVXR::sum_squares(A.star - Z.star %*% x))
    constraints <- list(CVXR::norm1(x[1:J]) <= QQ, x[1:J] >= lb)
    prob        <- CVXR::Problem(objective, constraints)
    sol         <- CVXR::solve(prob)
    
    b <- sol$getValue(x)
    names(b) <- colnames(Z)
    alert <- sol$status != "optimal"
    
    if (alert == TRUE) {
      stop(paste0("Estimation algorithm not converged! The algorithm returned the value:", 
                  sol$status, ". To check to what errors it corresponds go to 
                 'https://cvxr.rbind.io/cvxr_examples/cvxr_gentle-intro/'."))
    }
    
  } else {    # Optimization of all other cases
    
    if (p == "no norm") {
      res <-   nloptr::nloptr(x0          = rep(0, (J+KM)),
                              eval_f      = obj.fun.est,
                              lb          = c(rep(lb,J), rep(-Inf,KM)),
                              ub          = c(rep(Inf,J), rep(Inf,KM)),
                              opts        = opt.list,
                              Z = Z, V = V, A = A, J = J, QQ = QQ, KM = KM, p = pp)
    } else {
      if (dire == "==") {
        res <-   nloptr::nloptr(x0          = rep(0, (J+KM)),
                                eval_f      = obj.fun.est,
                                lb          = c(rep(lb,J), rep(-Inf,KM)),
                                ub          = c(rep(Inf,J), rep(Inf,KM)),
                                eval_g_eq   = norm.co.est, 
                                opts        = opt.list,
                                Z = Z, V = V, A = A, J = J, QQ = QQ, KM = KM, p = pp)
        
      } else if (dire == "<=") {
        res <-   nloptr::nloptr(x0          = rep(0, (J+KM)),
                                eval_f      = obj.fun.est,
                                lb          = c(rep(lb,J), rep(-Inf,KM)),
                                ub          = c(rep(Inf,J), rep(Inf,KM)),
                                eval_g_ineq = norm.co.est,
                                opts        = opt.list,
                                Z = Z, V = V, A = A, J = J, QQ = QQ, KM = KM, p = pp)
        
      }
    }  

    b         <- res$solution
    alert     <- res$status < 0 | res$status >= 5
  
  
    if (alert == TRUE) {
      stop(paste0("Estimation algorithm not converged! The algorithm returned the value:", 
                 res$status, ". To check to what errors it corresponds go to 
                 'https://nlopt.readthedocs.io/en/latest/NLopt_Reference/#return-values'."))
    }
  }
  
  names(b)  <- colnames(Z)
  
  return(b)
}


u.des.prep <- function(B, C, u.order, u.lags, coig.data, T0.tot, M, constant,
                       index, index.w, features, u.design, res) {
  
  ## Construct the polynomial terms in B
  if (u.order == 0) {                # Simple mean
    
    u.des.0 <- as.matrix(rep(1, T0.tot))
    
  } else if (u.order > 0) {         # Include covariates (u.order = 1 just covariates)
    
    if (coig.data == TRUE) {        # Take first differences of B and active covariates
      
      B.diff   <- NULL
      
      # Extract feature id from rownames of B
      feature.id <- unlist(purrr::map(stringr::str_split(rownames(B), "\\."), 1))
      
      # Create first differences feature-by-feature of the matrix B (not of C!!)
      for (feature in features) {
        BB       <- B[feature.id == feature,]
        B.diff   <- rbind(B.diff, BB - dplyr::lag(BB))
      }
      
      u.des.0 <- cbind(B.diff, C)[, index, drop=FALSE]    # combine with C
      
    } else if (coig.data == FALSE) {
      
      u.des.0 <- cbind(B, C)[, index, drop = FALSE]  # take active covariates
      
    }
    
    # Augment H with powers and interactions of B (not of C!!!)
    if (u.order > 1) {
      act.B <- sum(index.w)
      u.des.0 <- cbind(poly(u.des.0[,(1:act.B), drop = FALSE], degree = u.order, raw = TRUE, simple = TRUE),
                       u.des.0[,-(1:act.B), drop = FALSE])
    }
    
    # Include the constant if a global constant is not present
    # In case a constant is already specified lm.fit and qfit will automatically remove
    # the collinear covariates!!
    if (constant == FALSE) {
      u.des.0 <- cbind(u.des.0, rep(1, nrow(u.des.0)))
    }
  }

  ## Construct lags of B
  
  if (u.lags > 0) {
    B.lag <- NULL
    
    if (coig.data == TRUE) {
      # Take first differences of B
      B.diff   <- NULL
      
      # Extract feature id from rownames of B
      feature.id <- unlist(purrr::map(stringr::str_split(rownames(B), "\\."), 1))
      
      # Create first differences feature-by-feature of the matrix B (not of C!!)
      for (feature in features) {
        BB       <- B[feature.id == feature,]
        B.diff   <- rbind(B.diff, BB - dplyr::lag(BB))
      }
    }
    
    for (ll in seq_len(u.lags)) {
      B.l <- NULL
      for (feature in features) {
        if (coig.data == FALSE) {
          B.l <- rbind(B.l, dplyr::lag(B[feature.id == feature, , drop = FALSE], n = ll)[, index.w, drop = FALSE])
        } else {
          B.l <- rbind(B.l, dplyr::lag(B.diff[feature.id == feature, , drop = FALSE], n = ll)[, index.w, drop = FALSE])
        }
      }
      B.lag <- cbind(B.lag, B.l)
    }
    
    u.des.0 <- cbind(u.des.0, B.lag)
  }
  
  # If user provided check compatibility of the matrix and overwrite what has been created
  if (is.null(u.design) == FALSE) {
    if (is.matrix(u.design) == FALSE) {
      stop("The object u.design should be a matrix!!")
    }
    
    if (nrow(u.design) != nrow(res)) {
      stop(paste("The matrix u.design has", nrow(u.design),"rows when", nrow(res),
                 "where expected!"))
    }
    u.des.0 <- u.design
  }
  
  return(list(u.des.0 = u.des.0, f.id = feature.id))
}



e.des.prep <- function(B, C, P, e.order, e.lags, res, sc.pred, Y.donors, out.feat, features,
                       J, M, index, index.w, coig.data, T0, T0.tot, T1, constant, e.design) {
  
  
  # If the outcome variable is not among the features we need to create the
  # proper vector of residuals. Further, we force the predictors to be 
  # the outcome variable of the donors
  
  if (out.feat == FALSE) {
    e.res    <- sc.pred$data$Y.pre - sc.pred$est.results$Y.pre.fit
    
    if (coig.data == TRUE) {
      e.des.0 <- apply(Y.donors, 2, function(x) x - dplyr::lag(x))[, index.w]
      
      P.first <- P[1, ] - Y.donors[T0[1], ]
      P.diff  <- rbind(P.first, apply(P, 2, diff))[, index, drop = FALSE]
      e.des.1 <- P.diff
    } else {
      e.des.0  <- Y.donors[, index.w]
      e.des.1  <- P[, index, drop = FALSE]
    }
    
  } else if (out.feat == TRUE) {    # outcome variable is among features
    e.res <- res[1:T0[1],]
    
    ## Construct the polynomial terms in B (e.des.0) and P (e.des.1)
    if (e.order == 0) {
      
      e.des.0 <- as.matrix(rep(1, T0[1]))
      e.des.1 <- as.matrix(rep(1, T1))
      
    } else if (e.order > 0) {
      
      if (coig.data == TRUE) {
        
        ## Take first differences of B
        B.diff   <- NULL
        feature.id <- unlist(purrr::map(stringr::str_split(rownames(B), "\\."), 1))
        
        # Create first differences of the first feature (outcome) of the matrix B (not of C!!)
        BB       <- B[feature.id == features[1],]
        B.diff   <- rbind(B.diff, BB - dplyr::lag(BB))
        e.des.0 <- cbind(B.diff, C)[, index, drop=FALSE]

        ## Take first differences of P
        
        # Remove last observation of first feature from first period of P
        P.first <- c((P[1, (1:J)] - B[feature.id == features[1], , drop = FALSE][T0[1], ]),
                     P[1 ,-(1:J), drop = FALSE])
        
        
        # Take differences of other periods
        Pdiff    <- apply(P[, (1:J)], 2, diff)
        P.diff   <- rbind(P.first, cbind(Pdiff, P[-1,-(1:J)]))[, index, drop = FALSE]
        e.des.1  <- P.diff
        
      } else if (coig.data == FALSE) {
        e.des.0 <- cbind(B, C)[, index, drop = FALSE]
        e.des.1 <- P[, index, drop = FALSE]
      }

      # Augment H with powers and interactions of B (not of C!!!)
      if (e.order > 1) {
        act.B <- sum(index.w)
        e.des.0 <- cbind(poly(e.des.0[,(1:act.B), drop = FALSE], degree = e.order, raw = TRUE, simple = TRUE),
                         e.des.0[,-(1:act.B), drop = FALSE])
        e.des.1 <- cbind(poly(e.des.1[,(1:act.B), drop = FALSE], degree = e.order, raw = TRUE, simple = TRUE),
                         e.des.1[,-(1:act.B), drop = FALSE])
      }
      
      # Include the constant if a global constant is not present
      # In case a constant is already specified lm.fit will automatically remove
      # the collinear covariates!!
      if (constant == FALSE) {
        e.des.0 <- cbind(e.des.0, rep(1, nrow(e.des.0)))
        e.des.1 <- cbind(e.des.1, rep(1, nrow(e.des.1)))
      }
    }
    
    if (e.lags > 0) {
      # Construct lags of B and P
      B.lag <- NULL
      P.lag <- NULL
      
      # Take first differences of B and P
      B.diff   <- NULL
      feature.id <- unlist(purrr::map(stringr::str_split(rownames(B), "\\."), 1))
      
      # Create first differences of the first feature (outcome) of the matrix B (not of C!!)
      BB       <- B[feature.id == features[1], , drop = FALSE]
      B.diff   <- rbind(B.diff, BB - dplyr::lag(BB))

      ## Create first differences of P
      
      # Attach some pre-treatment value in order to avoid having missing values
      if (coig.data == FALSE) {
        PP <- rbind(B[feature.id == features[1], , drop = FALSE][((T0[1]-e.lags + 1):T0[1]), , drop = FALSE], 
                    P[ , (1:J), drop = FALSE])
      } else {
        PP <- rbind(B[feature.id == features[1], , drop = FALSE][((T0[1]-e.lags):T0[1]), , drop = FALSE], 
                    P[ , (1:J), drop = FALSE])
      }
      PP.diff <- PP - dplyr::lag(PP)

      for (ll in seq_len(e.lags)) {
        if (coig.data == FALSE) {
          P.l <- dplyr::lag(PP, n = ll)[, index.w, drop = FALSE][((e.lags+1):nrow(PP)), ,drop = FALSE]
        } else {
          P.l <- dplyr::lag(PP.diff, n = ll)[, index.w, drop = FALSE][((e.lags+2):nrow(PP)), ,drop = FALSE]
        }
        
        if (coig.data == FALSE) {
          B.l <- dplyr::lag(B[feature.id == features[1], , drop = FALSE], n = ll)[, index.w, drop = FALSE]
        } else {
          B.l <- dplyr::lag(B.diff[ , , drop = FALSE], n = ll)[, index.w, drop = FALSE]
        }
        
        B.lag <- cbind(B.lag, B.l)
        P.lag <- cbind(P.lag, P.l)
      }

      e.des.0 <- cbind(e.des.0, B.lag)
      e.des.1 <- cbind(e.des.1, P.lag)
    }
  }
  
  if (is.null(e.design) == FALSE) {
    if (is.matrix(e.design) == FALSE) {
      stop("The object e.design should be a matrix!!")
    }
    
    if (nrow(e.design) != nrow(e.res)) {
      stop(paste("The matrix e.design has", nrow(e.design),"rows when", nrow(e.res),
                 "where expected!"))
    }
    e.des.0 <- e.design
  }
  
  return(list(e.res = e.res, e.des.0 = e.des.0, e.des.1 = e.des.1))
}

DUflexGet <- function(u.des.0.na, C, f.id.na, M) {
  sel <- colnames(u.des.0.na) %in% colnames(C)
  D.b <- u.des.0.na[, !sel, drop = FALSE]
  D.c <- u.des.0.na[, sel, drop = FALSE]
  f.df <- data.frame(f.id.na)
  f.D <- fastDummies::dummy_cols(f.df, select_columns = "f.id", remove_selected_columns = TRUE)
  D.b.int <- matrix(NA, nrow = nrow(D.b), ncol=0)
  for (m in seq_len(M)) {
    D.b.int <- cbind(D.b.int, D.b*f.D[,m])
  }
  D <- cbind(D.b.int, D.c)
  
  return(D)
}


scpi.in <- function(xt, beta, Q, G, J, KM, p.int, QQ, dire, p, lb, w.lb.est, w.ub.est, opt.list) {
  
  # optimizing options
  opt.list <- prepareOptions(opt.list, p, dire, lb, "scpi")
  use.CVXR <- useCVXR(list(Q = QQ, p = p, dir = dire, lb = lb))

  # define optimization; min
  if (w.lb.est == TRUE) {
    
    if (use.CVXR == TRUE) { # handle L1 norm + inequality constraint
      a <- -2*G - 2*c(t(beta) %*% Q)
      d <- 2*sum(G*beta) + sum(beta*(Q %*% beta))
      
      x <- CVXR::Variable(J+KM)

      objective   <- CVXR::Minimize(-sum(CVXR::multiply(xt,x-beta)))
      
      if (lb == 0) {
        constraints <- list(CVXR::norm1(x[1:J]) <= QQ, 
                            x[1:J] >= 0,
                            CVXR::quad_form(x, Q) + sum(CVXR::multiply(a, x)) + d <= 0)
      } else {
        constraints <- list(CVXR::norm1(x[1:J]) <= QQ, 
                            CVXR::quad_form(x, Q) + sum(CVXR::multiply(a, x)) + d <= 0)
      }
      
      prob      <- CVXR::Problem(objective, constraints)
      prob_data <- CVXR::get_problem_data(prob, solver = "ECOS")
      ECOS_dims <- ECOS.dims_to_solver_dict(prob_data$data[["dims"]])
      solver_output <- ECOSolveR::ECOS_csolve(c = prob_data$data[["c"]],
                                              G = prob_data$data[["G"]],
                                              h = prob_data$data[["h"]],
                                              dims = ECOS_dims,
                                              A = prob_data$data[["A"]],
                                              b = prob_data$data[["b"]])
      sol      <- unpack_results(prob, solver_output, prob_data$chain, prob_data$inverse_data)
      alert    <- sol$status != "optimal"
      
      if (alert == TRUE) {
        lb.est <- NA
      } else {
        lb.est <- sol$value
      }

    } else {
      if (dire == "<=") {
        res.lb <-   nloptr(x0          = rep(0,(J+KM)),
                           eval_f      = obj.fun.min,
                           lb          = c(rep(lb,J), rep(-Inf,KM)),
                           ub          = c(rep(Inf,J), rep(Inf,KM)),
                           eval_g_ineq = double.ineq,
                           opts        = opt.list,
                           xt = xt, beta = beta, Q = Q, G = G, J = J, KM = KM, QQ = QQ, p.int = p.int)
        
        
      } else if (dire == "==") {
        res.lb <-   nloptr(x0          = rep(0,(J+KM)),
                           eval_f      = obj.fun.min,
                           lb          = c(rep(lb,J), rep(-Inf,KM)),
                           ub          = c(rep(Inf,J), rep(Inf,KM)),
                           eval_g_eq   = norm.equal,
                           eval_g_ineq = single.ineq,
                           opts        = opt.list,
                           xt = xt, beta = beta, Q = Q, G = G, J = J, KM = KM, QQ = QQ, p.int = p.int)
        
      } else if (dire == "NULL") {
        res.lb <-   nloptr(x0          = rep(0,(J+KM)),
                           eval_f      = obj.fun.min,
                           lb          = c(rep(lb,J), rep(-Inf,KM)),
                           ub          = c(rep(Inf,J), rep(Inf,KM)),
                           eval_g_ineq = single.ineq,
                           opts        = opt.list,
                           xt = xt, beta = beta, Q = Q, G = G, J = J, KM = KM, QQ = QQ, p.int = p.int)      
      }
  
      alert <- res.lb$status < 0 | res.lb$status >= 5
      flag  <- checkConstraints(res.lb, dire, 1.0e-2, 1.0e-2)  # allow for a little bit of slackness 
    
    
      if ((alert == TRUE) | (flag == TRUE)) {
        lb.est <- NA
      } else {
        lb.est <- res.lb$objective
      }
    }
  } else {
    lb.est <- NA
  }

  # define optimization; max
  if (w.ub.est == TRUE) {
    if (use.CVXR == TRUE) { # handle L1 norm + inequality constraint
      
      a <- -2*G - 2*c(t(beta) %*% Q)
      d <- 2*sum(G*beta) + sum(beta*(Q %*% beta))
      
      x <- CVXR::Variable(J+KM)
      
      objective   <- CVXR::Minimize(sum(CVXR::multiply(xt,x-beta)))
      
      if (lb == 0) {
        constraints <- list(CVXR::norm1(x[1:J]) <= QQ, 
                            x[1:J] >= 0,
                            CVXR::quad_form(x, Q) + sum(CVXR::multiply(a, x)) + d <= 0)
      } else {
        constraints <- list(CVXR::norm1(x[1:J]) <= QQ, 
                            CVXR::quad_form(x, Q) + sum(CVXR::multiply(a, x)) + d <= 0)
      }
      
      prob        <- CVXR::Problem(objective, constraints)
      prob        <- CVXR::Problem(objective, constraints)
      prob_data <- CVXR::get_problem_data(prob, solver = "ECOS")
      ECOS_dims <- ECOS.dims_to_solver_dict(prob_data$data[["dims"]])
      solver_output <- ECOSolveR::ECOS_csolve(c = prob_data$data[["c"]],
                                              G = prob_data$data[["G"]],
                                              h = prob_data$data[["h"]],
                                              dims = ECOS_dims,
                                              A = prob_data$data[["A"]],
                                              b = prob_data$data[["b"]])
      sol      <- unpack_results(prob, solver_output, prob_data$chain, prob_data$inverse_data)
      alert    <- sol$status != "optimal"

      if (alert == TRUE) {
        ub.est <- NA
      } else {
        ub.est <- -sol$value
      }
      
    } else {    
    
      if (dire == "<=") {
        res.ub <-   nloptr(x0          = rep(0,(J+KM)),
                           eval_f      = obj.fun.max,
                           lb          = c(rep(lb,J), rep(-Inf,KM)),
                           ub          = c(rep(Inf,J), rep(Inf,KM)),
                           eval_g_ineq = double.ineq,
                           opts        = opt.list,
                           xt = xt, beta = beta, Q = Q, G = G, J = J, KM = KM, QQ = QQ, p.int = p.int)
        
      } else if (dire == "==") {
        res.ub <-   nloptr(x0          = rep(0,(J+KM)),
                           eval_f      = obj.fun.max,
                           lb          = c(rep(lb,J), rep(-Inf,KM)),
                           ub          = c(rep(Inf,J), rep(Inf,KM)),
                           eval_g_eq   = norm.equal,
                           eval_g_ineq = single.ineq,
                           opts        = opt.list,
                           xt = xt, beta = beta, Q = Q, G = G, J = J, KM = KM, QQ = QQ, p.int = p.int)
        
      } else if (dire == "NULL") {
        res.ub <-   nloptr(x0          = rep(0,(J+KM)),
                           eval_f      = obj.fun.max,
                           lb          = c(rep(lb,J), rep(-Inf,KM)),
                           ub          = c(rep(Inf,J), rep(Inf,KM)),
                           eval_g_ineq = single.ineq,
                           opts        = opt.list,
                           xt = xt, beta = beta, Q = Q, G = G, J = J, KM = KM, QQ = QQ, p.int = p.int)      
      }
      
      alert <- res.ub$status < 0 | res.ub$status >= 5
      flag  <- checkConstraints(res.ub, dire, 1.0e-2, 1.0e-2)  # allow for a little bit of slackness 
      
      if ((alert == TRUE) | (flag == TRUE)) {
        ub.est <- NA
      } else {
        ub.est <- -res.ub$objective
      }
    }    
    
  } else {
    ub.est <- NA
  }
  
  return(c(lb.est, ub.est))
}

# prepare algorithm options
prepareOptions <- function(opt.list, p, dire, lb, input) {
  
  if (input == "scest") {
    if (is.null(opt.list$algorithm)) {
      if ((p == "L1") & (lb == -Inf)) {
        opt.list$algorithm <- 'NLOPT_LD_MMA'
      } else {
        opt.list$algorithm <- 'NLOPT_LD_SLSQP'
      } 
    }
    
    if (is.null(opt.list$xtol_rel)) opt.list$xtol_rel <- 1.0e-8
    if (is.null(opt.list$xtol_abs)) opt.list$xtol_abs <- 1.0e-8
    if (is.null(opt.list$ftol_rel)) opt.list$ftol_rel <- 1.0e-8
    if (is.null(opt.list$ftol_abs)) opt.list$ftol_abs <- 1.0e-8
    if (is.null(opt.list$maxeval))  opt.list$maxeval  <- 5000
    
    if (dire == "==") {
      if (is.null(opt.list$tol_constraints_eq)) {
        opt.list$tol_constraints_eq <- 1.0e-8
      } else {
        opt.list$tol_constraints_eq <- opt.list$tol_constraints_eq[1] 
      }
    } else if (dire == "<=") {
      if (is.null(opt.list$tol_constraints_ineq)) {
        opt.list$tol_constraints_ineq <- 1.0e-8
      } else {
        opt.list$tol_constraints_ineq <- opt.list$tol_constraints_ineq[1] 
      }
    }
  }

  if (input == "scpi") {
    if (is.null(opt.list$algorithm)) {
        opt.list$algorithm <- 'NLOPT_LD_SLSQP'
      } 
    
    if (is.null(opt.list$xtol_rel)) opt.list$xtol_rel <- 1.0e-8
    if (is.null(opt.list$xtol_abs)) opt.list$xtol_abs <- 1.0e-8
    if (is.null(opt.list$ftol_rel)) opt.list$ftol_rel <- 1.0e-8
    if (is.null(opt.list$ftol_abs)) opt.list$ftol_abs <- 1.0e-8
    if (is.null(opt.list$maxeval))  opt.list$maxeval  <- 5000
    
    if (dire == "NULL") {
      if (is.null(opt.list$tol_constraints_ineq)) {
        opt.list$tol_constraints_ineq <- 1.0e-8
      } else {
        opt.list$tol_constraints_ineq <- opt.list$tol_constraints_ineq[1] 
      }
    } else if (dire == "==") {
      if (is.null(opt.list$tol_constraints_ineq)) {
        opt.list$tol_constraints_ineq <- 1.0e-8
      } else {
        opt.list$tol_constraints_ineq <- opt.list$tol_constraints_ineq[1] 
      }
      if (is.null(opt.list$tol_constraints_eq)) {
        opt.list$tol_constraints_eq <- 1.0e-8
      } else {
        opt.list$tol_constraints_eq <- opt.list$tol_constraints_eq[1] 
      }
    } else if (dire == "<=") {
      if (is.null(opt.list$tol_constraints_ineq)) {
        opt.list$tol_constraints_ineq <- rep(1.0e-8,2)
      } else {
        opt.list$tol_constraints_ineq <- rep(opt.list$tol_constraints_ineq[1], 2) 
      }
    }
  }
  return(opt.list) 
}


# function to check that inequality and equality constraints are satisfied
checkConstraints <- function(nloptr.obj, dir, tol_eq, tol_ineq) {
  
  if (dir == "NULL") {
    flag <- nloptr.obj$eval_g_ineq(nloptr.obj$solution)$constraints > tol_ineq
    
  } else if (dir == "==") {
    flag1 <- nloptr.obj$eval_g_ineq(nloptr.obj$solution)$constraints    > tol_ineq
    flag2 <- abs(nloptr.obj$eval_g_eq(nloptr.obj$solution)$constraints) > tol_eq
    
    flag <- flag1 | flag2
    
  } else if(dir == "<=") {
    flag1 <- nloptr.obj$eval_g_ineq(nloptr.obj$solution)$constraints[1] > tol_ineq
    flag2 <- nloptr.obj$eval_g_ineq(nloptr.obj$solution)$constraints[2] > tol_ineq
    
    flag <- flag1 | flag2
  }
  return(flag)
}

# Prediction interval, for e
scpi.out <- function(res, x, eval, e.method, alpha, e.lb.est, e.ub.est, verbose) {

  neval <- nrow(eval)
  e.1 <- e.2 <- lb <- ub <- NA

  if (e.lb.est == TRUE | e.ub.est == TRUE) {

    if (e.method == "gaussian") {
      x.more   <- rbind(eval, x)
      fit      <- predict(y=res, x=x, eval=x.more, type="lm", verbose = verbose)
      e.mean   <- fit[1:neval]
      res.fit  <- fit[-(1:neval)]
      var.pred <- predict(y=log((res-res.fit)^2), x=x, eval=x.more, type="lm", verbose = verbose)
      e.sig2   <- exp(var.pred[1:neval])

      eps <- sqrt(-log(alpha)*2*e.sig2)
      lb <- e.mean - eps
      ub <- e.mean + eps

      # save mean and variance of u, only for sensitivity analysis
      e.1 <- e.mean
      e.2 <- e.sig2

    } else if (e.method == "ls") {
      x.more  <- rbind(eval, x)
      fit     <- predict(y=res, x=x, eval=x.more, type="lm", verbose = verbose)
      e.mean  <- fit[1:neval]
      res.fit <- fit[-(1:neval)]

      var.pred <- predict(y=log((res-res.fit)^2), x=x, eval=x.more, type="lm", verbose = verbose)
      e.sig    <- sqrt(exp(var.pred[1:neval]))
      res.st   <- (res-res.fit)/sqrt(exp(var.pred[-(1:neval)]))

      lb <- e.mean + e.sig * quantile(res.st, alpha)
      ub <- e.mean + e.sig * quantile(res.st, 1-alpha)
      
      # save mean and variance of u, only for sensitivity analysis
      e.1 <- e.mean
      e.2 <- e.sig^2
      
    } else if (e.method == "qreg") {
      e.pred  <- predict(y=res, x=x, eval=eval, type="qreg", tau=c(alpha, 1-alpha), verbose = verbose)
      lb <- e.pred[,1]
      ub <- e.pred[,2]

    }
  }

  return(list(lb = lb, ub = ub, e.1 = e.1, e.2 = e.2))
}



# square root matrix
sqrtm <- function(A) {
  decomp <- svd(A)
  decomp$d[decomp$d < 0] <- 0
  rootA  <- decomp$u %*% diag(sqrt(decomp$d)) %*% t(decomp$u)
  return(rootA)
}

# conditional prediction
predict <- function(y, x, eval, type="lm", tau=NULL, verbose) {
  
  if ((nrow(x) <= ncol(x)) & verbose) {
    warning("Consider specifying a less complicated model for e. The number of observations used
         to parametrically predict moments is smaller than the number of covariates used. Consider reducing either the number
         of lags (e.lags) or the order of the polynomial (e.order)!")
  }
  
  if (type == "lm") {
    betahat <- .lm.fit(x, y)$coeff

  } else if (type == "qreg") {
    if (is.null(tau)) {
      tau <- c(0.05, 0.95)
    }

    betahat <- rrq(y~x-1, tau=tau)$coeff
  }
  pred <- eval %*% betahat
  return(pred)
}


## Auxiliary function that estimates degrees of freedom

df.EST <- function(w.constr, w, A, B, J, KM){
  if ((w.constr[["name"]] == "ols") | (w.constr[["p"]] == "no norm")) {
    df <- J 
    
  } else if ((w.constr[["name"]] == "lasso") | ((w.constr[["p"]] == "L1") & (w.constr[["dir"]] == "<="))) {
    df <- sum(w != 0) 
    
  } else if ((w.constr[["name"]] == "simplex") | ((w.constr[["p"]] == "L1") & (w.constr[["dir"]] == "=="))) {
    df <- sum(w != 0) - 1 
  
  } else if ((w.constr[["name"]] == "ridge") | (w.constr[["p"]] == "L2")) {
    d <- svd(B)$d
    d[d < 0] <- 0
    df <- sum(d^2/(d^2+w.constr[["lambda"]]))

  } 
    
  # add degrees of freedom coming from C block
  df <- df + KM
  
  return(df)
}





u.sigma.est <- function(u.mean, u.sigma, res, Z, V, index, TT, M, df) {
  
  if      (u.sigma == "HC0") { # White (1980)
    vc <- 1
  }
  
  else if (u.sigma == "HC1") { # MacKinnon and White (1985)
    vc <- TT/(TT-df)
  }
  
  else if (u.sigma == "HC2") { # MacKinnon and White (1985)
    PP <- Z %*% base::solve(t(Z)%*% V %*% Z) %*% t(Z) %*% V
    vc <- 1/(1-diag(PP))
  }
  
  else if (u.sigma == "HC3") { # Davidson and MacKinnon (1993)
    PP <- Z %*% base::solve(t(Z)%*% V %*% Z) %*% t(Z) %*% V
    vc <- 1/(1-diag(PP))^2
  }
  
  else if (u.sigma == "HC4") { # Cribari-Neto (2004)
    PP <- Z %*% base::solve(t(Z)%*% V %*% Z) %*% t(Z) %*% V
    CN <- as.matrix((TT*M)*diag(PP)/df)
    dd <- apply(CN, 1, function(x) min(4,x))
    vc <- as.matrix(NA, length(res), 1)
    for (ii in seq_len(length(res))) {
      vc[ii] <- 1/(1 - diag(PP)[ii])^dd[ii]
    }
  }
  
  Omega  <- diag(c((res-u.mean)^2)*vc)
  Sigma  <- t(Z) %*% V %*% Omega %*% V %*% Z / (TT^2)
  
  return(list(Omega = Omega, Sigma = Sigma))
}


local.geom <- function(w.constr, rho, rho.max, res, B, C, coig.data, T0.tot, J, w, verbose) {
  
  Q   <- w.constr[["Q"]]
  if (is.character(rho)) rho <- regularize.w(rho, rho.max, res, B, C, coig.data, T0.tot)
  
  if ((w.constr[["name"]] == "simplex") | ((w.constr[["p"]] == "L1") & (w.constr[["dir"]] == "=="))) {
    index.w <- abs(w) > rho
    index.w <- regularize.check(w, index.w, rho, verbose)
    
    w.star  <- w
    w.star[!index.w] <- 0
    
    Q.star  <- sum(w.star)  
    

  } else if ((w.constr[["name"]] == "lasso") | ((w.constr[["p"]] == "L1") & (w.constr[["dir"]] == "<="))) {
    
    if ((sum(abs(w)) >= Q - rho*sqrt(J)) & (sum(abs(w)) <= Q)) {
      Q.star <- sum(abs(w))
    } else {
      Q.star <- Q
    }
    index.w <- abs(w) > rho
    index.w <- regularize.check(w, index.w, rho, verbose)
    
    w.star  <- w
    w.star[!index.w] <- 0    
    
  } else if( (w.constr[["name"]] == "ridge") | (w.constr[["p"]] == "L2")){

    if (sqrt((sum(w^2)) >= Q - rho) & (sqrt(sum(w^2)) <= Q)) {
      Q.star <- sqrt(sum(w^2))
    } else {
      Q.star <- Q
    }
    
    index.w <- rep(TRUE, length(w))
    w.star  <- w
    
  } else {
    Q.star  <- Q
    w.star  <- w
    index.w <- rep(TRUE, length(w))
  }

  w.constr[["Q"]] <- Q.star
  return(list(w.constr = w.constr, w.star = w.star, index.w = index.w, rho = rho, Q.star = Q.star)) 
}


regularize.w <- function(rho, rho.max, res, B, C, coig.data, T0.tot) {
  
  if (rho == "type-1") {
    sigma.u  <- sqrt(mean((res-mean(res))^2))
    sigma.bj <- min(apply(B, 2, sd))
    CC       <- sigma.u/sigma.bj
    
  } else if (rho == "type-2"){
    sigma.u   <- sqrt(mean((res-mean(res))^2))
    sigma.bj2 <- min(apply(B, 2, var))
    sigma.bj  <- max(apply(B, 2, sd))
    CC        <- sigma.bj*sigma.u/sigma.bj2
    
  } else if (rho == "type-3"){
    sigma.bj2 <- min(apply(B, 2, var))
    sigma.bju <- max(apply(B, 2, function(bj) cov(bj, res)))
    CC        <- sigma.bju/sigma.bj2
  }
  
    
  if (coig.data == TRUE) { # cointegration
    c <- 1
  } else {        # iid or ar
    c <- 0.5
  }
  
  rho <- (CC*(log(T0.tot))^c)/(sqrt(T0.tot))
  
  if (is.null(rho.max) == FALSE)  rho <- min(rho, rho.max)         
  
  return(rho)
}

regularize.check <- function(w, index.w, rho, verbose) {
  if (sum(index.w) == 0) {
    index.w <- rank(-w) <= 1
    if (verbose){
      cat("Warning: regularization paramater was too high (", round(rho, digits = 3), "). ", sep = "")
      cat("We set it so that at least one component in w is non-zero.")
    }
  }
  return(index.w)
}


useCVXR <- function(w.constr) {
  flag <- FALSE
  if ((w.constr$p == "L1") & (w.constr$dir == "<=")) flag <- TRUE
  return(flag)
}

executionTime <- function(T0, J, T1, sims, cores, name){
  Ttot <- sum(T0)
  tincr <- Ttot/1000
  coefsJ <- c(-0.54755616,0.09985644)
  
  time <- sum(c(1,J)*coefsJ)    # Time for J donors, T0 = 1000, sims = 10
  time <- time*sims/10          # rescale for simulations
  time <- time/cores            # rescale by number of cores
  time <- time*tincr            # rescale for number of obs
  time <- time*T1               # rescale by post intervention periods
  
  time <- time/60               # Convert in minutes
  time <- ceiling(time)         # Take smallest larger integer
  time <- time*2
  
  if (name == "lasso") {
    time <- time*8
  }
  
  if (time < 1) {
    toprint <- "Maximum expected execution time: less than a minute.\n"
  } else if (time == 1) {
    toprint <- paste0("Maximum expected execution time: ",time," minute.\n")
  } else {
    toprint <- paste0("Maximum expected execution time: ",time," minutes.\n")
  }
  cat(toprint)
  cat("\n")
}



