###############################################################################
### Auxiliary functions for estimation

# Quadratic loss function
obj.fun.est <- function(x, Z, V, A, J, QQ, KM, p) {
  f <- x %*% t(Z) %*% V %*% Z %*% x - 2 * t(A) %*% V %*% Z %*% x
  g <- 2 * t(Z) %*% V %*% Z %*% x - 2 * t(t(A) %*% V %*% Z)

  return(list("objective" = f,
              "gradient"  = g))
}

obj.fun.est.sr <- function(x, Z, V, A, J, Q1, Q2, KMI, S) {
  f <- x %*% t(Z) %*% V %*% Z %*% x - 2 * t(A) %*% V %*% Z %*% x
  g <- 2 * t(Z) %*% V %*% Z %*% x - 2 * t(t(A) %*% V %*% Z)
  
  return(list("objective" = f,
              "gradient"  = g))
}

obj.fun.est.multi <- function(x, Z, V, A, J, QQ, KMI, p, S) {
  f <- x %*% t(Z) %*% V %*% Z %*% x - 2 * t(A) %*% V %*% Z %*% x
  g <- 2 * t(Z) %*% V %*% Z %*% x - 2 * t(t(A) %*% V %*% Z)

  return(list("objective" = f,
              "gradient"  = g))
}


## Constraint on the norm

# Single treated unit
norm.co.est <- function(x, Z, V, A, J, QQ, KM, p) {

  if (p == 1) {
    av <- rep(1, J)
    av[x[1:J] < 0] <- -1
    ja <- c(av, rep(0, KM))
    co <- sum(abs(x[1:J])) - QQ
  } else {
    ja <- c(2 * x[1:J], rep(0, KM))
    co <- sum(x[1:J]^2) - QQ^2
  }

  return(list("constraints" = co,
              "jacobian"    = ja))
}

# Functions for L1-L2

norm.L1 <- function(x, Z, V, A, J, Q1, Q2, KMI, S) {
  
  av <- rep(1, J)
  av[x[1:J] < 0] <- -1
  ja <- c(av, rep(0, KMI))
  co <- S %*% abs(x[1:J]) - Q1
  
  ja.vec <- matrix(ja, nrow = nrow(S), ncol = length(ja), byrow = TRUE)
  ja.vec[, 1:J] <- ja.vec[, 1:J] * S
  
  return(list("constraints" = co,
              "jacobian"    = ja.vec))  
}

norm.L2 <- function(x, Z, V, A, J, Q1, Q2, KMI, S) {
  
  ja <- c(2 * x[1:J], rep(0, KMI))
  co <- S %*% x[1:J]^2 - Q2^2
  
  ja.vec <- matrix(ja, nrow = nrow(S), ncol = length(ja), byrow = TRUE)
  ja.vec[, 1:J] <- ja.vec[, 1:J] * S
  
  return(list("constraints" = co,
              "jacobian"    = ja.vec))
}

# Multiple treated units
norm.co.est.multi <- function(x, Z, V, A, J, QQ, KMI, p, S) {

  if (p == 1) {
    av <- rep(1, J)
    av[x[1:J] < 0] <- -1
    ja <- c(av, rep(0, KMI))
    co <- S %*% abs(x[1:J]) - QQ
  } else {
    ja <- c(2 * x[1:J], rep(0, KMI))
    co <- S %*% x[1:J]^2 - QQ^2
  }

  ja.vec <- matrix(ja, nrow = nrow(S), ncol = length(ja), byrow = TRUE)
  ja.vec[, 1:J] <- ja.vec[, 1:J] * S
  
  return(list("constraints" = co,
              "jacobian"    = ja.vec))
}


### Auxiliary functions for inference

# Prepare objective functions
obj.fun.min <- function(x, xt, beta, Q, G, J, KMI, QQ, p.int, S) {
  f <- -sum(xt * (x - beta))
  g <- -xt

  return(list("objective" = f,
              "gradient"  = g))
}

obj.fun.max <- function(x, xt, beta, Q, G, J, KMI, QQ, p.int, S) {
  f <- sum(xt * (x - beta))
  g <- xt

  return(list("objective" = f,
              "gradient"  = g))
}

obj.fun.min.sr <- function(x, xt, beta, Q, G, J, KMI, Q1, Q2, S) {
  f <- -sum(xt * (x - beta))
  g <- -xt
  
  return(list("objective" = f,
              "gradient"  = g))
}

obj.fun.max.sr <- function(x, xt, beta, Q, G, J, KMI, Q1, Q2, S) {
  f <- sum(xt * (x - beta))
  g <- xt
  
  return(list("objective" = f,
              "gradient"  = g))
}


# Prepare inequality constraint(s): loss functions constraint + (inequality norm constraint)

# Unique inequality constraint
single.ineq <- function(x, xt, beta, Q, G, J, KMI, QQ, p.int, S) {
  a <- -2 * G - 2 * c(t(beta) %*% Q)
  d <- 2 * sum(G * beta) + sum(beta * (Q %*% beta))

  co <- x %*% Q %*% x + sum(a * x) + d
  ja <- 2 * Q %*% x + a

  return(list("constraints" = co,
              "jacobian"    = ja))
}

# Inequality constraints: loss function + L1-L2 norm
double.ineq <- function(x, xt, beta, Q, G, J, KMI, QQ, p.int, S) {
  # Loss function constraint
  a <- -2 * G - 2 * c(t(beta) %*% Q)
  d <- 2 * sum(G * beta) + sum(beta * (Q %*% beta))

  co1 <- x %*% Q %*% x + sum(a * x) + d
  ja1 <- 2 * Q %*% x + a

  # Norm constraint
  if (p.int == 1) {
    co2 <- S %*% abs(x[1:J]) - QQ
    av <- rep(1, J)
    av[x[1:J] < 0] <- -1
    ja <- c(av, rep(0, KMI))
  } else {
    co2 <- S %*% (x[1:J]^2) - QQ^2
    ja <- c(2 * x[1:J], rep(0, KMI))
  }

  ja.vec <- matrix(ja, nrow = nrow(S), ncol = length(ja), byrow = TRUE)
  ja.vec[, 1:J] <- ja.vec[, 1:J] * S
  ja <- c(rbind(t(ja1), ja.vec))    # vectorize matrix row-by-row

  return(list("constraints" = c(co1,co2),
              "jacobian"    = ja))
}

# Eventual equality constraint on norm
norm.equal <- function(x, xt, beta, Q, G, J, KMI, QQ, p.int, S) {

  if (p.int == 1) {
    co <- S %*% abs(x[1:J]) - QQ
    av <- rep(1, J)
    av[x[1:J] < 0] <- -1
    ja <- c(av, rep(0, KMI))
  } else {
    co <- S %*% (x[1:J]^2) - QQ^2
    ja <- c(2 * x[1:J], rep(0, KMI))
  }

  ja.vec <- matrix(ja, nrow = nrow(S), ncol = length(ja), byrow = TRUE)
  ja.vec[, 1:J] <- ja.vec[, 1:J] * S

  return(list("constraints" = co,
              "jacobian"    = ja.vec))
}

norm.equal.sr <- function(x, xt, beta, Q, G, J, KMI, Q1, Q2, S) {
  
  co <- S %*% abs(x[1:J]) - Q1
  av <- rep(1, J)
  av[x[1:J] < 0] <- -1
  ja <- c(av, rep(0, KMI))

  ja.vec <- matrix(ja, nrow = nrow(S), ncol = length(ja), byrow = TRUE)
  ja.vec[, 1:J] <- ja.vec[, 1:J] * S
  
  return(list("constraints" = co,
              "jacobian"    = ja.vec))
}

double.ineq.sr <- function(x, xt, beta, Q, G, J, KMI, Q1, Q2, S) {
  # Loss function constraint
  a <- -2 * G - 2 * c(t(beta) %*% Q)
  d <- 2 * sum(G * beta) + sum(beta * (Q %*% beta))
  
  co1 <- x %*% Q %*% x + sum(a * x) + d
  ja1 <- 2 * Q %*% x + a
  
  co2 <- S %*% (x[1:J]^2) - Q2^2
  ja <- c(2 * x[1:J], rep(0, KMI))

  ja.vec <- matrix(ja, nrow = nrow(S), ncol = length(ja), byrow = TRUE)
  ja.vec[, 1:J] <- ja.vec[, 1:J] * S
  ja <- c(rbind(t(ja1), ja.vec))    # vectorize matrix row-by-row
  
  return(list("constraints" = c(co1,co2),
              "jacobian"    = ja))
}


# Auxiliary function that creates the constraints to be passed to the optimization problem
w.constr.OBJ <- function(w.constr, A, Z, V, J, KM, M) {
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

      feature.id <- unlist(purrr::map(stringr::str_split(rownames(Z), "\\."), 2))

      Qfeat <- c()
      for (feat in unique(feature.id)) {
        Af <- A[feature.id == feat, , drop=FALSE]
        Zf <- Z[feature.id == feat, , drop=FALSE]
        Vf <- V[feature.id == feat, feature.id == feat, drop=FALSE]

        if (nrow(Af) >= 5) {
          QQ <- tryCatch({
                            aux <- shrinkage.EST("ridge", Af, Zf, Vf, J, KM)
                            Q <- aux$Q
                          }, warning={}, error={}, finally={})
          Qfeat <- c(Qfeat, QQ)
        }
      }

      if (is.null(Qfeat)) Qfeat <- shrinkage.EST("ridge", A, Z, V, J, KM)$Q
      w.constr[["Q"]]      <- min(Qfeat, na.rm=TRUE)
      w.constr[["lambda"]] <- aux$lambda
    }

    w.constr <- list(lb     = -Inf,
                     p      = "L2",
                     dir    = "<=",
                     Q      = w.constr[["Q"]],
                     name   = "ridge",
                     lambda = w.constr[["lambda"]])

  } else if (w.constr[["name"]] == "L1-L2") {

    if (!("Q2" %in% names(w.constr))) {
      
      feature.id <- unlist(purrr::map(stringr::str_split(rownames(Z), "\\."), 2))
      Qfeat <- c()
      for (feat in unique(feature.id)) {
        Af <- A[feature.id == feat, , drop=FALSE]
        Zf <- Z[feature.id == feat, , drop=FALSE]
        Vf <- V[feature.id == feat, feature.id == feat, drop=FALSE]
        
        if (nrow(Af) >= 5) {
          QQ <- tryCatch({
            aux <- shrinkage.EST("ridge", Af, Zf, Vf, J, KM)
            Q2 <- aux$Q
          }, warning={}, error={}, finally={})
          Qfeat <- c(Qfeat, QQ)
        }
      }
      
      if (is.null(Qfeat)) Qfeat <- shrinkage.EST("ridge", A, Z, V, J, KM)$Q
      w.constr[["Q2"]]      <- max(Qfeat, na.rm=TRUE)
      w.constr[["lambda"]] <- aux$lambda
    }
  
    w.constr <- list(lb     = -Inf,
                     p      = "L1-L2",
                     dir    = "==/<=",
                     Q      = 1,
                     Q2     = w.constr[["Q2"]],
                     name   = "L1-L2",
                     lambda = w.constr[["lambda"]])
    
  } else {
    # if constraint is entirely user specified just check everything is fine
    if (!(all(c('p','dir','Q','lb') %in% names(w.constr)))) {
      stop("If 'name' is not specified, w.constr should be a list whose elements 
            must be named 'p','dir','Q','lb'.")
    }

    if (!(w.constr[["p"]] %in% c("no norm","L1","L2","L1-L2"))) {
      stop("Specify either p = 'no norm' (no constraint on the norm of weights),
          p = 'L1' (L1-norm), p = 'L2' (L2-norm)")
    } else if (w.constr[["p"]] == "no norm") {
      w.constr[["dir"]] <- "NULL"
    }

    if (!(w.constr[["dir"]] %in% c("<=","==","==/<=","NULL"))) {
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
    lambd   <- sig.wls^2*(J + KM) / sum(wls$coef^2, na.rm = TRUE)           # rule of thumb for lambda (Hoerl et al, 1975)
    Q       <- sqrt(sum(wls$coef^2, na.rm = TRUE)) / (1 + lambd)            # convert lambda into Q
    
    if (is.nan(Q) | (nrow(Z) <= ncol(Z) + 10)) { # reduce dimensionality of the problem if more params than obs
      lasso.cols <- b.est(A, Z, J, KM, list(dir = "<=", lb = -Inf, p = "L1", Q = 1), V, NULL)
      active.cols <- abs(lasso.cols) > 1e-8
      if (sum(active.cols) >= (max(nrow(A) - 10, 2)) ) {
        active.cols <-  rank(-abs(lasso.cols)) <= max(nrow(A) - 10, 2)
      }
      Z.sel <- Z[,active.cols,drop=FALSE]
      wls     <- lm(A ~ Z.sel - 1, weights = diag(V))
      sig.wls <- sigma(wls)
      lambd   <- sig.wls^2*(ncol(Z.sel)+KM)/sum(wls$coef^2, na.rm = TRUE)     # rule of thumb for lambda (Hoerl et al, 1975)
      Q       <- sqrt(sum(wls$coef^2, na.rm = TRUE)) / (1 + lambd)            # convert lambda into Q
    }
  }

  return(list(Q = Q, lambda = lambd))
}


# Auxiliary function that solves the (un)constrained problem to estimate b
# depending on the desired method
b.est <- function(A, Z, J, KM, w.constr, V, opt.list) {

  dire <- w.constr[["dir"]]
  lb   <- w.constr[["lb"]]
  p    <- w.constr[["p"]]
  QQ   <- w.constr[["Q"]]

  if (p == "no norm") pp <- 0
  if (p == "L1") pp <- 1
  if (p == "L2") pp <- 2
  if (p == "L1-L2") {
    pp <- NULL
    Q2 <- w.constr[["Q2"]]
  }    

  opt.list <- prepareOptions(opt.list, p, dire, lb, "scest")
  use.CVXR <- useCVXR(w.constr)

  if (use.CVXR == TRUE) { # handle L1 norm + inequality constraint

    x <- CVXR::Variable(J + KM)

    objective   <- CVXR::Minimize(CVXR::quad_form(A - Z %*% x, V))
    constraints <- list(CVXR::norm1(x[1:J]) <= QQ, x[1:J] >= lb)
    prob        <- CVXR::Problem(objective, constraints)
    sol         <- CVXR::solve(prob)

    b <- sol$getValue(x)
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
      
    } else if (p == "L1-L2") {
      S <- matrix(1, nrow = 1, ncol = J)
      res <-   nloptr::nloptr(x0          = rep(0, (J+KM)),
                              eval_f      = obj.fun.est.sr,
                              lb          = c(rep(lb,J), rep(-Inf,KM)),
                              ub          = c(rep(Inf,J), rep(Inf,KM)),
                              eval_g_eq   = norm.L1, 
                              eval_g_ineq = norm.L2, 
                              opts        = opt.list,
                              Z = Z, V = V, A = A, J = J, Q1 = QQ, Q2 = Q2, KMI = KM, S = S)
      
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

  if (is.matrix(b)) {
    b <- b[, 1, drop = TRUE]
  }
  names(b) <- colnames(Z)

  return(b)
}

# Auxiliary function that solves the (un)constrained problem to estimate b
# depending on the desired method - Multiple treated units case
b.est.multi <- function(A, Z, J, KMI, I, w.constr, V, opt.list, S) {

  # The constraint is symmetric in the shape across treated units (J, KM, Q might change)
  dire  <- w.constr[[1]]$dir
  lb    <- w.constr[[1]]$lb
  p     <- w.constr[[1]]$p
  QQ    <- unlist(lapply(w.constr, function(x) x$Q))

  if (p == "no norm") pp <- 0
  if (p == "L1") pp <- 1
  if (p == "L1-L2") {
    pp <- NULL
    Q2 <- unlist(lapply(w.constr, function(x) x$Q2))
  }  
  
  opt.list <- prepareOptions(opt.list, p, dire, lb, "scest", I)
  use.CVXR <- useCVXR(w.constr[[1]])
  
  Jtot <- sum(unlist(J))
  
  if (use.CVXR == TRUE) { # handle L1 norm + inequality constraint
    x <- CVXR::Variable(Jtot+KMI)
    
    objective   <- CVXR::Minimize(CVXR::quad_form(A - Z %*% x, V))
    constraints <- list(x[1:Jtot] >= lb)
    j.lb <- 1
    for (i in seq_len(I)) {
      j.ub <- j.lb + J[[i]] - 1 
      constraints <- append(constraints, list(CVXR::norm1(x[j.lb:j.ub]) <= QQ[i]))
      j.lb <- j.ub + 1
    }
    prob        <- CVXR::Problem(objective, constraints)
    sol         <- CVXR::solve(prob)
    
    b <- sol$getValue(x)
    alert <- sol$status != "optimal"
    
    if (alert == TRUE) {
      stop(paste0("Estimation algorithm not converged! The algorithm returned the value:", 
                  sol$status, ". To check to what errors it corresponds go to 
                  'https://cvxr.rbind.io/cvxr_examples/cvxr_gentle-intro/'."))
    }
    
  } else {    # Optimization of all other cases
    
    if (p == "no norm") {
      res <-   nloptr::nloptr(x0          = rep(0, (Jtot+KMI)),
                              eval_f      = obj.fun.est.multi,
                              lb          = c(rep(lb,Jtot), rep(-Inf,KMI)),
                              ub          = c(rep(Inf,Jtot), rep(Inf,KMI)),
                              opts        = opt.list,
                              Z = Z, V = V, A = A, J = Jtot, QQ = QQ, KMI = KMI, p = pp, S = S)
      
    } else if (p == "L1-L2") {
      res <-   nloptr::nloptr(x0          = rep(0, (Jtot+KMI)),
                              eval_f      = obj.fun.est.sr,
                              lb          = c(rep(lb,Jtot), rep(-Inf,KMI)),
                              ub          = c(rep(Inf,Jtot), rep(Inf,KMI)),
                              eval_g_eq   = norm.L1, 
                              eval_g_ineq = norm.L2, 
                              opts        = opt.list,
                              Z = Z, V = V, A = A, J = Jtot, Q1 = QQ, Q2 = Q2, KMI = KMI, S = S)
      
    } else {
      if (dire == "==") {
        res <-   nloptr::nloptr(x0          = rep(0, (Jtot+KMI)),
                                eval_f      = obj.fun.est.multi,
                                lb          = c(rep(lb,Jtot), rep(-Inf,KMI)),
                                ub          = c(rep(Inf,Jtot), rep(Inf,KMI)),
                                eval_g_eq   = norm.co.est.multi, 
                                opts        = opt.list,
                                Z = Z, V = V, A = A, J = Jtot, QQ = QQ, KMI = KMI, p = pp, S = S)
        
      } else if (dire == "<=") {
        res <-   nloptr::nloptr(x0          = rep(0, (Jtot+KMI)),
                                eval_f      = obj.fun.est.multi,
                                lb          = c(rep(lb,Jtot), rep(-Inf,KMI)),
                                ub          = c(rep(Inf,Jtot), rep(Inf,KMI)),
                                eval_g_ineq = norm.co.est.multi,
                                opts        = opt.list,
                                Z = Z, V = V, A = A, J = Jtot, QQ = QQ, KMI = KMI, p = pp, S = S)
        
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

  if (is.matrix(b)) {
    rownames(b)  <- colnames(Z)
  } else {
    names(b) <- colnames(Z)
  }

  return(b)
}


V.prep <- function(type, B, T0.features, I) {
  if (type == "separate") { # Default (separate fit)
    V <- diag(dim(B)[1])
    
  } else if (type == "pooled") {
    
    dim.V <- unlist(lapply(T0.features, function(x) sum(unlist(x))))
    max.dim <- max(dim.V)
    ones <- matrix(1, nrow = I, ncol = 1)
    eye <- diag(1, nrow = max.dim, ncol = max.dim)
    V <- kronecker(ones %*% t(ones), eye) # structure if T0*M was balanced across treated unit
    
    sel <- matrix(TRUE, nrow = nrow(V), ncol = ncol(V))
    for (i in seq_len(I)) { # trim V according to length of pre-treatment period
      if (dim.V[i] < max.dim) {
        shift <- (i - 1) * max.dim
        sel[(shift + dim.V[i] + 1) : (shift + max.dim), ] <- FALSE
        sel[, (shift + dim.V[i] + 1) : (shift + max.dim)] <- FALSE
      }
    }
    row.trim <- rowSums(sel) != 0
    col.trim <- colSums(sel) != 0
    
    V <- V[row.trim, col.trim] / I^2
  }

  rownames(V) <- rownames(B)
  colnames(V) <- rownames(B)

  return(V)
}



u.des.prep <- function(B, C, u.order, u.lags, coig.data, T0.tot, constant,
                       index, index.w, features, feature.id, u.design, res, verbose) {

  ## Construct the polynomial terms in B
  if (u.order == 0) {                # Simple mean

    u.des.0 <- as.matrix(rep(1, T0.tot))

  } else if (u.order > 0) {         # Include covariates (u.order = 1 just covariates)

    if (coig.data == TRUE) {        # Take first differences of B and active covariates

      B.diff   <- NULL

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
      name.tr <- lapply(strsplit(rownames(u.des.0), "\\."), "[[", 1)[[1]]
      act.B <- sum(index.w)
      u.des.poly <- poly(u.des.0[, (1:act.B), drop = FALSE], degree = u.order, raw = TRUE, simple = TRUE)
      colnames(u.des.poly) <- paste(name.tr, colnames(u.des.poly), sep = ".")
      u.des.0 <- cbind(u.des.poly,
                       u.des.0[, -(1:act.B), drop = FALSE])

    }

    # Include the constant if a global constant is not present
    # In case a constant is already specified lm.fit and qfit will automatically remove
    # the collinear covariates!!
    if (constant == FALSE) {
      u.des.0 <- cbind(u.des.0, rep(1, nrow(u.des.0)))
      name.tr <- lapply(strsplit(rownames(u.des.0), "\\."), "[[", 1)[[1]]
      colnames(u.des.0) <- c(colnames(u.des.0[, -ncol(u.des.0), drop = FALSE]),
                             paste0(name.tr, ".0.constant"))
    }
  }

  ## Construct lags of B
  if (u.lags > 0) {

    B.lag <- NULL
    if (coig.data == TRUE) {
      # Take first differences of B
      B.diff   <- NULL

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
          B.l <- rbind(B.l, dplyr::lag(B[feature.id == feature, , drop = FALSE], n = ll))
        } else {
          B.l <- rbind(B.l, dplyr::lag(B.diff[feature.id == feature, , drop = FALSE], n = ll))
        }
      }
      B.lag <- cbind(B.lag, B.l)
    }
    name.tr <- lapply(strsplit(rownames(u.des.0), "\\."), "[[", 1)[[1]]
    colnames(B.lag) <- rep(paste0(name.tr,".lag"), ncol(B.lag))
    u.des.0 <- cbind(u.des.0, B.lag[, index.w, drop = FALSE])
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

  return(list(u.des.0 = u.des.0))
}


e.des.prep <- function(B, C, P, e.order, e.lags, res, sc.pred, Y.donors, out.feat, features,
                       J, index, index.w, coig.data, T0, T1, constant, e.design, P.diff.pre) {

  # If the outcome variable is not among the features we need to create the
  # proper vector of residuals. Further, we force the predictors to be
  # the outcome variable of the donors
  sel <- c()

  aux <- trendRemove(P)
  C <- trendRemove(C)$mat
  index <- index[aux$sel]
  P <- aux$mat
  if (!is.null(P.diff.pre)) P.diff.pre <- trendRemove(as.matrix(P.diff.pre))$mat
  
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
    e.res <- res[1:T0[1], , drop = FALSE]

    ## Construct the polynomial terms in B (e.des.0) and P (e.des.1)
    if (e.order == 0) {

      e.des.0 <- as.matrix(rep(1, T0[1]))
      e.des.1 <- as.matrix(rep(1, T1))

    } else if (e.order > 0) {
      feature.id <- unlist(purrr::map(stringr::str_split(rownames(B), "\\."), 2))

      if (coig.data == TRUE) {

        ## Take first differences of B
        B.diff   <- NULL

        # Create first differences of the first feature (outcome) of the matrix B (not of C!!)
        BB       <- B[feature.id == features[1], ]
        B.diff   <- rbind(B.diff, BB - dplyr::lag(BB))
        e.des.0 <- cbind(B.diff, C[feature.id == features[1], ])[, index, drop = FALSE]


        ## Take first differences of P
        # Remove last observation of first feature from first period of P
        P.first <- c((P[1, (1:J), drop = FALSE] - B[feature.id == features[1], , drop = FALSE][T0[1], ]),
                     P[1 , -(1:J), drop = FALSE])


        # Take differences of other periods
        if (nrow(P) > 2) {
          Pdiff    <- apply(P[, (1:J), drop = FALSE], 2, diff)
          P.diff   <- rbind(P.first, cbind(Pdiff, P[-1, -(1:J), drop = FALSE]))[, index, drop = FALSE]
        } else if (nrow(P) == 2) {
          Pdiff    <- t(as.matrix(apply(P[, (1:J), drop = FALSE], 2, diff)))
          P.diff   <- rbind(P.first, cbind(Pdiff, P[-1, -(1:J), drop = FALSE]))[, index, drop = FALSE]
        } else {
          P.diff <- matrix(P.first, 1, length(P.first))[, index, drop = FALSE]
        }
        e.des.1  <- P.diff

      } else {
        e.des.0 <- cbind(B, C)[feature.id == features[1], index, drop = FALSE]
        e.des.1 <- P[, index, drop = FALSE]
      }
      
      # Augment H with powers and interactions of B (not of C!!!)
      if (e.order > 1) {
        act.B <- sum(index.w)
        e.des.0 <- cbind(poly(e.des.0[,(1:act.B), drop = FALSE], degree = e.order, raw = TRUE, simple = TRUE),
                         e.des.0[, -(1:act.B), drop = FALSE])
        e.des.1 <- cbind(poly(e.des.1[,(1:act.B), drop = FALSE], degree = e.order, raw = TRUE, simple = TRUE),
                         e.des.1[, -(1:act.B), drop = FALSE])
      }

      # Include the constant if a global constant is not present
      # In case a constant is already specified lm.fit will automatically remove
      # the collinear covariates!!
      if (constant == FALSE) {
        e.des.0 <- cbind(e.des.0, rep(1, nrow(e.des.0)))
        e.des.1 <- cbind(e.des.1, rep(1, nrow(e.des.1)))
      }
    }
  
    nolag <- FALSE
    if (is.null(P.diff.pre) == FALSE) {
      e.des.1 <- P.diff.pre[, index, drop = FALSE]
      nolag <- TRUE
    }

    if (e.lags > 0 && nolag == FALSE) {
      # Construct lags of B and P
      B.lag <- NULL
      P.lag <- NULL

      # Take first differences of B and P
      B.diff   <- NULL
      feature.id <- unlist(purrr::map(stringr::str_split(rownames(B), "\\."), 2))

      # Create first differences of the first feature (outcome) of the matrix B (not of C!!)
      BB       <- B[feature.id == features[1], , drop = FALSE]
      B.diff   <- rbind(B.diff, BB - dplyr::lag(BB))

      ## Create first differences of P
      # Attach some pre-treatment value in order to avoid having missing values
      if (coig.data == FALSE) {
        PP <- rbind(B[feature.id == features[1], , drop = FALSE][((T0[1] - e.lags + 1):T0[1]), , drop = FALSE],
                    P[, (1:J), drop = FALSE])
      } else {
        PP <- rbind(B[feature.id == features[1], , drop = FALSE][((T0[1] - e.lags):T0[1]), , drop = FALSE],
                    P[, (1:J), drop = FALSE])
      }
      PP.diff <- PP - dplyr::lag(PP)

      for (ll in seq_len(e.lags)) {
        if (coig.data == FALSE) {
          P.l <- dplyr::lag(PP, n = ll)[, index.w, drop = FALSE][((e.lags+1):nrow(PP)), , drop = FALSE]
        } else {
          P.l <- dplyr::lag(PP.diff, n = ll)[, index.w, drop = FALSE][((e.lags+2):nrow(PP)), , drop = FALSE]
        }

        if (coig.data == FALSE) {
          B.l <- dplyr::lag(B[feature.id == features[1], , drop = FALSE], n = ll)[, index.w, drop = FALSE]
        } else {
          B.l <- dplyr::lag(B.diff[, , drop = FALSE], n = ll)[, index.w, drop = FALSE]
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
      stop(paste("The matrix e.design has", nrow(e.design), "rows when", nrow(e.res),
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
  D.b.int <- matrix(NA, nrow = nrow(D.b), ncol = 0)
  for (m in seq_len(ncol(f.D))) {
    D.b.int <- cbind(D.b.int, D.b*f.D[,m])
  }

  D <- cbind(D.b.int, D.c)

  return(D)
}


insampleUncertaintyGet <- function(Z.na, V.na, P.na, beta, S, Sigma.root, J, KMI, I,
                                   w.constr.inf, Q.star, Q2.star, lb, TT, sims, cores, verbose, 
                                   w.lb.est, w.ub.est, opt.list.inf) {

  Q <- t(Z.na) %*% V.na %*% Z.na / TT
  colnames(Q) <- colnames(Z.na)

  if (w.constr.inf[["p"]] == "no norm") p.int <- 0
  if (w.constr.inf[["p"]] == "L1") p.int <- 1
  if (w.constr.inf[["p"]] == "L2") p.int <- 2
  if (w.constr.inf[["p"]] == "L1-L2") p.int <- NULL
  
  # Algorithm initial value is lower bound unless -Inf
  x0 <- lb
  x0[is.infinite(lb)] <- 0

  jj <- nrow(P.na)

  # optimizing options
  
  opt.list <- prepareOptions(opt.list.inf, w.constr.inf[["p"]], w.constr.inf[["dir"]], lb, "scpi", I)
  use.CVXR <- useCVXR(list(Q = Q.star, p = w.constr.inf[["p"]], dir = w.constr.inf[["dir"]], lb = lb))
  Jtot <- sum(unlist(J))

  iters <- round(sims / 10)
  perc  <- 0
  # simulate

  if (cores == 1) {
    vsig <- matrix(NA, nrow = sims, ncol = 2 * jj)

    for (sim in seq_len(sims)) {
      rem <- sim %% iters
      if ((rem == 0) && verbose) {
        perc <- perc + 10
        cat(paste(sim, "/", sims, " iterations completed (", perc, "%)", " \r", sep = ""))
        utils::flush.console()
      }

      zeta    <- rnorm(length(beta))
      G       <- Sigma.root %*% zeta

      for (hor in seq_len(jj)) {
        xt <- P.na[hor, ]

        output  <- scpi.in(xt = xt, beta = beta, Q = Q, G = G, J = J, KMI = KMI, I = I, Jtot = Jtot, S = S,
                           p.int = p.int, QQ = Q.star, QQ2 = Q2.star, lb = lb, x0 = x0,
                           dire = w.constr.inf[["dir"]], p = w.constr.inf[["p"]],
                           w.lb.est = w.lb.est, w.ub.est = w.ub.est, opt.list = opt.list, use.CVXR = use.CVXR)

        vsig[sim, hor]      <- output[1]
        vsig[sim, hor + jj] <- output[2]

      }
    }

  } else if (cores >= 1) {

    progress <- function(n) {
      rem <- n %% iters
      if ((rem == 0) && verbose) {
        perc <- n/sims * 100
        cat(paste(n, "/", sims, " iterations completed (", perc, "%)", " \r", sep = ""))
        utils::flush.console()
      }
    }
    opts <- list(progress=progress)

    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)

    vsig <- foreach::foreach(i = 1 : sims,
                             .packages = c('nloptr','CVXR'),
                             .export   = c('scpi.in','obj.fun.min','obj.fun.max', 'checkConstraints', 'useCVXR',
                                           'single.ineq','double.ineq', 'norm.equal', 'prepareOptions',
                                           'obj.fun.min.sr','obj.fun.max.sr','double.ineq.sr', 'norm.equal.sr'),
                             .combine  = rbind,
                             .options.snow = opts) %dorng% {
                               
                               zeta   <- rnorm(length(beta))
                               G      <- Sigma.root %*% zeta
                               
                               ub.sim <- c()
                               lb.sim <- c()
                               
                               for (hor in seq_len(jj)) {
                                 xt <- P.na[hor,]
                                 
                                 output  <- scpi.in(xt = xt, beta = beta, Q = Q, G = G, J = J, KMI = KMI, I = I, Jtot = Jtot, S = S,
                                                    p.int = p.int, QQ = Q.star, QQ2 = Q2.star, lb = lb, x0 = x0,
                                                    dire = w.constr.inf[["dir"]], p = w.constr.inf[["p"]], 
                                                    w.lb.est = w.lb.est, w.ub.est = w.ub.est, opt.list = opt.list, use.CVXR = use.CVXR)
                                 
                                 lb.sim      <- append(lb.sim, output[1])
                                 ub.sim      <- append(ub.sim, output[2])
                                 
                               }
                               
                               c(lb.sim, ub.sim)
                             }
    
    parallel::stopCluster(cl)
  }
  
  return(vsig)
}


scpi.in <- function(xt, beta, Q, G, J, KMI, I, Jtot, S, p.int, QQ, QQ2, dire, p, lb, x0, 
                    w.lb.est, w.ub.est, opt.list, use.CVXR) {
  # define optimization; min
  if (w.lb.est == TRUE) {
    
    if (use.CVXR == TRUE) { # handle L1 norm + inequality constraint
      a <- -2 * G - 2 * c(t(beta) %*% Q)
      d <- 2 * sum(G * beta) + sum(beta * (Q %*% beta))

      x <- CVXR::Variable(Jtot+KMI)

      objective   <- CVXR::Minimize(-sum(CVXR::multiply(xt,x - beta)))
      
      constraints <- list(CVXR::quad_form(x, Q) + sum(CVXR::multiply(a, x)) + d <= 0)
      
      if (lb[1] > - Inf) {
        constraints <- append(constraints, list(x[1:Jtot] >= lb))
      }
      
      j.lb <- 1
      for (i in seq_len(I)) {
        j.ub <- j.lb + J[[i]] - 1 
        constraints <- append(constraints, list(CVXR::norm1(x[j.lb:j.ub]) <= QQ[i]))
        j.lb <- j.ub + 1
      }      
      
      prob     <- CVXR::Problem(objective, constraints)
      sol      <- CVXR::solve(prob)
      alert    <- !(sol$status %in% c("optimal","optimal_inaccurate"))
      
      if (alert == TRUE) {
        lb.est <- NA
      } else {
        lb.est <- sol$value
      }
      
    } else {

      if (dire == "<=") {
        res.lb <-   nloptr(x0          = c(x0, rep(0,KMI)),
                           eval_f      = obj.fun.min,
                           lb          = c(lb, rep(-Inf,KMI)),
                           ub          = c(rep(Inf,Jtot), rep(Inf,KMI)),
                           eval_g_ineq = double.ineq,
                           opts        = opt.list,
                           xt = xt, beta = beta, Q = Q, G = G, J = Jtot, KMI = KMI, 
                           QQ = QQ, p.int = p.int, S = S)
  
      } else if (dire == "==") {
        res.lb <-   nloptr(x0          = c(x0, rep(0,KMI)),
                           eval_f      = obj.fun.min,
                           lb          = c(lb, rep(-Inf,KMI)),
                           ub          = c(rep(Inf,Jtot), rep(Inf,KMI)),
                           eval_g_eq   = norm.equal,
                           eval_g_ineq = single.ineq,
                           opts        = opt.list,
                           xt = xt, beta = beta, Q = Q, G = G, J = Jtot, KMI = KMI, 
                           QQ = QQ, p.int = p.int, S = S)

      } else if (dire == "==/<=") {
        res.lb <-   nloptr(x0          = c(x0, rep(0,KMI)),
                           eval_f      = obj.fun.min.sr,
                           lb          = c(lb, rep(-Inf,KMI)),
                           ub          = c(rep(Inf,Jtot), rep(Inf,KMI)),
                           eval_g_eq   = norm.equal.sr,
                           eval_g_ineq = double.ineq.sr,
                           opts        = opt.list,
                           xt = xt, beta = beta, Q = Q, G = G, J = Jtot, KMI = KMI, 
                           Q1 = QQ, Q2 = QQ2, S = S)

      } else if (dire == "NULL") {
        res.lb <-   nloptr(x0          = c(x0, rep(0,KMI)),
                           eval_f      = obj.fun.min,
                           lb          = c(lb, rep(-Inf,KMI)),
                           ub          = c(rep(Inf,Jtot), rep(Inf,KMI)),
                           eval_g_ineq = single.ineq,
                           opts        = opt.list,
                           xt = xt, beta = beta, Q = Q, G = G, J = Jtot, KMI = KMI, 
                           QQ = QQ, p.int = p.int, S = S)      
      }
      
      alert <- res.lb$status < 0 | res.lb$status >= 5
      #flag  <- checkConstraints(res.lb, dire, 1.0e-2, 1.0e-2)  # allow for a little bit of slackness 
      
      
      if ((alert == TRUE) ) { #| (flag == TRUE)
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
      
      x <- CVXR::Variable(Jtot+KMI)
      
      objective   <- CVXR::Minimize(sum(CVXR::multiply(xt,x-beta)))
      
      constraints <- list(CVXR::quad_form(x, Q) + sum(CVXR::multiply(a, x)) + d <= 0)
      
      if (!is.infinite(lb[1])) {
        constraints <- append(constraints, list(x[1:J] >= 0))
      }
      
      j.lb <- 1
      for (i in seq_len(I)) {
        j.ub <- j.lb + J[[i]] - 1 
        constraints <- append(constraints, list(CVXR::norm1(x[j.lb:j.ub]) <= QQ[i]))
        j.lb <- j.ub + 1
      }      
      
      prob     <- CVXR::Problem(objective, constraints)
      sol      <- CVXR::solve(prob)
      alert    <- !(sol$status %in% c("optimal","optimal_inaccurate"))
      
      if (alert == TRUE) {
        ub.est <- NA
      } else {
        ub.est <- -sol$value
      }
      
    } else {    
      
      if (dire == "<=") {
        res.ub <-   nloptr(x0          = c(x0, rep(0,KMI)),
                           eval_f      = obj.fun.max,
                           lb          = c(lb, rep(-Inf,KMI)),
                           ub          = c(rep(Inf,Jtot), rep(Inf,KMI)),
                           eval_g_ineq = double.ineq,
                           opts        = opt.list,
                           xt = xt, beta = beta, Q = Q, G = G, J = Jtot, KMI = KMI, 
                           QQ = QQ, p.int = p.int, S = S)
        
      } else if (dire == "==") {
        res.ub <-   nloptr(x0          = c(x0, rep(0,KMI)),
                           eval_f      = obj.fun.max,
                           lb          = c(lb, rep(-Inf,KMI)),
                           ub          = c(rep(Inf,Jtot), rep(Inf,KMI)),
                           eval_g_eq   = norm.equal,
                           eval_g_ineq = single.ineq,
                           opts        = opt.list,
                           xt = xt, beta = beta, Q = Q, G = G, J = Jtot, KMI = KMI, 
                           QQ = QQ, p.int = p.int, S = S)

      } else if (dire == "==/<=") {
        res.ub <-   nloptr(x0          = c(x0, rep(0,KMI)),
                           eval_f      = obj.fun.max.sr,
                           lb          = c(lb, rep(-Inf,KMI)),
                           ub          = c(rep(Inf,Jtot), rep(Inf,KMI)),
                           eval_g_eq   = norm.equal.sr,
                           eval_g_ineq = double.ineq.sr,
                           opts        = opt.list,
                           xt = xt, beta = beta, Q = Q, G = G, J = Jtot, KMI = KMI, 
                           Q1 = QQ, Q2 = QQ2, S = S)        
        
      } else if (dire == "NULL") {
        res.ub <-   nloptr(x0          = c(x0, rep(0,KMI)),
                           eval_f      = obj.fun.max,
                           lb          = c(lb, rep(-Inf,KMI)),
                           ub          = c(rep(Inf,Jtot), rep(Inf,KMI)),
                           eval_g_ineq = single.ineq,
                           opts        = opt.list,
                           xt = xt, beta = beta, Q = Q, G = G, J = Jtot, KMI = KMI, 
                           QQ = QQ, p.int = p.int, S = S)      
      }
      
      alert <- res.ub$status < 0 | res.ub$status >= 5
      #flag  <- checkConstraints(res.ub, dire, 1.0e-2, 1.0e-2)  # allow for a little bit of slackness 
      
      if ((alert == TRUE)) { #| (flag == TRUE)
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
prepareOptions <- function(opt.list, p, dire, lb, input, I = 1) {
  
  if (input == "scest") {
    if (is.null(opt.list$algorithm)) {
      if ((p == "L1") && (any(lb == -Inf))) {
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
        opt.list$tol_constraints_eq <- rep(1.0e-8, I)
      } else {
        opt.list$tol_constraints_eq <- rep(opt.list$tol_constraints_eq[1], I) 
      }
      
    } else if (dire == "<=") {
      if (is.null(opt.list$tol_constraints_ineq)) {
        opt.list$tol_constraints_ineq <- rep(1.0e-8, I)
      } else {
        opt.list$tol_constraints_ineq <- rep(opt.list$tol_constraints_ineq[1], I) 
      }
      
    } else if (dire == "==/<=") {
      opt.list$ftol_rel <- 1.0e-32
      opt.list$ftol_abs <- 1.0e-32
      if (is.null(opt.list$tol_constraints_eq)) {
        opt.list$tol_constraints_eq <- rep(1.0e-8, I)
      } else {
        opt.list$tol_constraints_eq <- rep(opt.list$tol_constraints_eq[1], I) 
      }
      if (is.null(opt.list$tol_constraints_ineq)) {
        opt.list$tol_constraints_ineq <- rep(1.0e-8, I)
      } else {
        opt.list$tol_constraints_ineq <- rep(opt.list$tol_constraints_ineq[1], I) 
      } 
    }
  }
  
  if (input == "scpi") {
    
    if (is.null(opt.list$algorithm)) {
      if ((p == "L1") && (any(lb == -Inf))) {
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
        opt.list$tol_constraints_eq <- rep(1.0e-8, I)
      } else {
        opt.list$tol_constraints_eq <- rep(opt.list$tol_constraints_eq[1], I) 
      }
      
    } else if (dire == "<=") {
      opt.list$ftol_rel <- 1.0e-32
      opt.list$ftol_abs <- 1.0e-32
      if (is.null(opt.list$tol_constraints_ineq)) {
        opt.list$tol_constraints_ineq <- rep(1.0e-8, I + 1)
      } else {
        opt.list$tol_constraints_ineq <- rep(opt.list$tol_constraints_ineq[1], I + 1) 
      }
      
    } else if (dire == "==/<=") {
      opt.list$ftol_rel <- 1.0e-32
      opt.list$ftol_abs <- 1.0e-32
      if (is.null(opt.list$tol_constraints_ineq)) {
        opt.list$tol_constraints_ineq <- rep(1.0e-8, I + 1)
      } else {
        opt.list$tol_constraints_ineq <- rep(opt.list$tol_constraints_ineq[1], I + 1)
      }
      if (is.null(opt.list$tol_constraints_eq)) {
        opt.list$tol_constraints_eq <- rep(1.0e-8, I)
      } else {
        opt.list$tol_constraints_eq <- rep(opt.list$tol_constraints_eq[1], I) 
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
    
    flag <- flag1 | any(flag2)
    
  } else if(dir == "<=") {
    flag <- nloptr.obj$eval_g_ineq(nloptr.obj$solution)$constraints > tol_ineq
    flag <- any(flag)
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
      fit      <- predict(y=res, x=x, eval=x.more, type="lm")
      e.mean   <- fit[1:neval]
      res.fit  <- fit[-(1:neval)]
      
      var.pred <- predict(y=log((res-res.fit)^2), x=x, eval=x.more, type="lm")
      e.sig2   <- exp(var.pred[1:neval])

      q.pred <- predict(y=res-res.fit, x=x, eval=x.more, type="qreg", tau = c(0.25,0.75))
      IQ.pred <- q.pred[1:neval,2] - q.pred[1:neval,1]
      IQ.pred <- abs(IQ.pred)
      e.sig <- apply(cbind(sqrt(e.sig2), IQ.pred/1.34), 1, min)
      eps <- sqrt(-log(alpha)*2)*e.sig
      
      lb <- e.mean - eps
      ub <- e.mean + eps
      
      # save mean and variance of u, only for sensitivity analysis
      e.1 <- e.mean
      e.2 <- e.sig^2
      
    } else if (e.method == "ls") {
      x.more  <- rbind(eval, x)
      fit     <- predict(y=res, x=x, eval=x.more, type="lm")
      e.mean  <- fit[1:neval]
      res.fit <- fit[-(1:neval)]
      
      var.pred <- predict(y=log((res-res.fit)^2), x=x, eval=x.more, type="lm")
      e.sig    <- sqrt(exp(var.pred[1:neval]))
      res.st   <- (res-res.fit)/sqrt(exp(var.pred[-(1:neval)]))
      
      q.pred <- predict(y=res-res.fit, x=x, eval=x.more, type="qreg", tau = c(0.25,0.75))
      IQ.pred <- q.pred[1:neval,2] - q.pred[1:neval,1]
      IQ.pred <- abs(IQ.pred)
      e.sig <- apply(cbind(e.sig, IQ.pred/1.34), 1, min)
      
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


simultaneousPredGet <- function(vsig, T1, T1.tot, I, u.alpha, e.alpha, e.res.na, e.des.0.na, e.des.1,
                                w.lb.est, w.ub.est, w.bounds, w.name) {

  vsigUB <- vsig[, (T1.tot + 1):(2 * T1.tot), drop = FALSE]
  vsigLB <- vsig[, 1:T1.tot, drop = FALSE]

  pi.e   <- scpi.out(res = e.res.na, x = e.des.0.na, eval = e.des.1, 
                     e.method = "gaussian", alpha = e.alpha/2, e.lb.est = TRUE, e.ub.est =  TRUE)

  w.lb.joint <- w.ub.joint <- c()
  
  j.min <- 1
  
  for (i in seq_len(I)) {
    j.max <- T1[[i]] + j.min - 1

#     if (w.name %in% c("ols", "ridge", "L1-L2")) {
#       lb.joint <- quantile(apply(vsigLB[, j.min:j.max, drop = FALSE], 1, quantile, na.rm = TRUE, probs = 0.01),
# 						   na.rm = TRUE, probs = u.alpha/2)
#       ub.joint <- quantile(apply(vsigUB[, j.min:j.max, drop = FALSE], 1, quantile, na.rm = TRUE, probs = 0.99),
# 						   na.rm = TRUE, probs = (1-u.alpha/2))
#       
#     } else {
    lb.joint <- quantile(apply(vsigLB[, j.min:j.max, drop = FALSE], 1, min, na.rm = TRUE), probs = u.alpha/2)
    ub.joint <- quantile(apply(vsigUB[, j.min:j.max, drop = FALSE], 1, max, na.rm = TRUE), probs = (1-u.alpha/2))
#    }
    
    w.lb.joint <- c(w.lb.joint, rep(lb.joint, T1[[i]]))
    w.ub.joint <- c(w.ub.joint, rep(ub.joint, T1[[i]]))
    j.min <- j.max + 1
    
  }

  eps <- 1
  if (length(pi.e$e.1) > 1) {
    eps <- c()
    for (i in seq_len(I)) {
      eps <- c(eps, rep(sqrt(log(T1[[i]] + 1)), T1[[i]]))
    }
  }
  
  e.lb.joint <- pi.e$lb * eps
  e.ub.joint <- pi.e$ub * eps
  
  if (w.lb.est == FALSE) w.lb.joint <- w.bounds[, 1]
  if (w.ub.est == FALSE) w.ub.joint <- w.bounds[, 2]
  
  MU <- e.ub.joint + w.ub.joint 
  ML <- e.lb.joint + w.lb.joint
  
  return(list(MU=MU, ML=ML))
}

epskappaGet <- function(P, rho.vec, beta, I, joint = FALSE) {
  P.list <- mat2list(P)
  beta.list <- mat2list(beta)
  
  epskappa <- c()
  for (i in seq_len(I)) {
    epskappai <- apply(P.list[[i]], 1, 
                       function(x) rho.vec[[i]]^2*sum(abs(x))/(2*sqrt(sum(beta.list[[i]]^2))))  
    epskappa <- c(epskappa, epskappai)
  }
  
  if (joint == TRUE) epskappa <- max(epskappa)
  
  return(epskappa)
}

# square root matrix
sqrtm <- function(A) {
  decomp <- svd(A)
  decomp$d[decomp$d < 0] <- 0
  rootA  <- decomp$u %*% diag(sqrt(decomp$d)) %*% t(decomp$u)
  return(rootA)
}

# conditional prediction
predict <- function(y, x, eval, type="lm", tau=NULL, verbose = FALSE) {
  
  if (type == "lm") {
    betahat <- .lm.fit(x, y)$coeff
    pred <- eval %*% betahat
    
  } else if (type == "qreg") {
    
    tryCatch(
      {
        betahat <- Qtools::rrq(y~x-1, tau=tau)$coefficients
      },
      
      warning = function(war) {
        message("Warning produced when estimating moments of the out-of-sample residuals with quantile regressions.")
        war$call <- NULL
        if (verbose) warning(war)
      }
    )
    
    betahat <- suppressWarnings(Qtools::rrq(y~x-1, tau=tau)$coefficients)
    pred <- eval %*% betahat
  }
  
  
  return(pred)
}


## Auxiliary function that estimates degrees of freedom

df.EST <- function(w.constr, w, B, J, KM){
  if ((w.constr[["name"]] == "ols") || (w.constr[["p"]] == "no norm")) {
    df <- J 
    
  } else if ((w.constr[["name"]] == "lasso") || ((w.constr[["p"]] == "L1") && (w.constr[["dir"]] == "<="))) {
    df <- sum(abs(w) >= 1e-6) 
    
  } else if ((w.constr[["name"]] == "simplex") || ((w.constr[["p"]] == "L1") && (w.constr[["dir"]] == "=="))) {
    df <- sum(abs(w) >= 1e-6) - 1
    
  } else if ((w.constr[["name"]] == "ridge") || (w.constr[["name"]] == "L1-L2") || (w.constr[["p"]] == "L2")) {
    d <- svd(B)$d
    d[d < 0] <- 0
    df <- sum(d^2/(d^2+w.constr[["lambda"]]))
    
  } 
  
  # add degrees of freedom coming from C block
  df <- df + KM
  
  return(df)
}


u.sigma.est <- function(u.mean, u.sigma, res, Z, V, index, TT, df) {
  
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
    CN <- as.matrix((TT)*diag(PP)/df)
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
  Q2.star <- NULL
  
  if (is.character(rho)) rho <- regularize.w(rho, rho.max, res, B, C, coig.data, T0.tot)
  
  if ((w.constr[["name"]] == "simplex") || ((w.constr[["p"]] == "L1") && (w.constr[["dir"]] == "=="))) {
    index.w <- abs(w) > rho
    index.w <- regularize.check(w, index.w, rho, verbose)
    w.star  <- w
    w.star[!index.w] <- 0
    
    Q.star  <- sum(w.star)  
    
  } else if ((w.constr[["name"]] == "lasso") || ((w.constr[["p"]] == "L1") && (w.constr[["dir"]] == "<="))) {
    
    if ((sum(abs(w)) >= Q - rho*sqrt(J)) && (sum(abs(w)) <= Q)) {
      Q.star <- sum(abs(w))
    } else {
      Q.star <- Q
    }
    index.w <- abs(w) > rho
    index.w <- regularize.check(w, index.w, rho, verbose)
    
    w.star  <- w
    w.star[!index.w] <- 0    
    
  } else if ((w.constr[["name"]] == "ridge") || (w.constr[["p"]] == "L2")) {
    
    if (sqrt((sum(w^2)) >= Q - rho) && (sqrt(sum(w^2)) <= Q)) {
      Q.star <- sqrt(sum(w^2))
    } else {
      Q.star <- Q
    }

    index.w <- rep(TRUE, length(w))
    w.star  <- w
    
  } else if (w.constr[["name"]] == "L1-L2") {
    
    index.w <- abs(w) > rho
    index.w <- regularize.check(w, index.w, rho, verbose)
    w.star  <- w
    w.star[!index.w] <- 0
    
    Q.star  <- sum(w.star)  

    if (sqrt((sum(w^2)) >= Q - rho) && (sqrt(sum(w^2)) <= Q)) {
      Q2.star <- sqrt(sum(w^2))
    } else {
      Q2.star <- w.constr[["Q2"]] 
    }
    w.constr[["Q2"]] <- Q2.star
    
  } else {
    Q.star  <- Q
    w.star  <- w
    index.w <- rep(TRUE, length(w))
  }
  
  w.constr[["Q"]] <- Q.star
  
  return(list(w.constr = w.constr, w.star = w.star, index.w = index.w, rho = rho, Q.star = Q.star, Q2.star = Q2.star)) 
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
      warning(paste0("Regularization paramater was too high (", round(rho, digits = 3), "). ",
                     "We set it so that at least one component in w is non-zero."), immediate. = TRUE, call. = FALSE)
    }
  }
  return(index.w)
}

local.geom.2step <- function(w, r, rho.vec, w.constr, Q, I) {
  beta <- c(w, r)
  w.list <- mat2list(as.matrix(w))
  
  ## Constraint on the norm of the weights
  if (w.constr[[1]]$p == "no norm") { # Unconstrained problem
    rhoj.vec <- rho.vec # auxiliary list never used in this case
    
  } else if (w.constr[[1]]$p == "L1") { 
    rhoj.vec <- rho.vec
    w.norm <- unlist(lapply(w.list, function(x) sum(abs(x))))
    
  } else if (w.constr[[1]]$p %in% c("L2","L1-L2")) {
    rhoj.vec <- c()
    for (i in seq_len(I)) {
      rhoj.vec[i] <- 2*sum(abs(w.list[[i]]))*rho.vec[i]
    }
    w.norm <- unlist(lapply(w.list, function(x) sum(x^2)))
  }
  
  # Check if constraint is equality or inequality
  if (w.constr[[1]]$dir %in% c("<=", "==/<=")) {
    active <- 1*((w.norm - Q) > -rhoj.vec)
    Q <- active*(w.norm - Q) + Q 
  }
  
  ## Constraint on lower bound of the weights
  lb <- c()
  for (i in seq_len(I)) {
    if (w.constr[[i]]$lb == 0) {
      active <- 1*(w.list[[i]] < rhoj.vec[[i]])
      lb <- c(lb, rep(0, length(w.list[[i]])) + active*w.list[[i]])
    } else {
      lb <- c(lb, rep(-Inf, length(w.list[[i]])))
    }
  }
  
  return(list(Q = Q, lb = lb))
}


useCVXR <- function(w.constr) {
  flag <- FALSE
  if ((w.constr$p == "L1") && (w.constr$dir == "<=")) flag <- TRUE
  return(flag)
}

executionTime <- function(T0, J, I, T1, sims, cores, name){
  Ttot <- sum(T0)
  tincr <- Ttot/1000
  coefsJ <- c(-0.54755616,0.09985644)
  
  time <- sum(c(1,J)*coefsJ)    # Time for J donors, T0 = 1000, sims = 10
  time <- ceiling(time)*sims/10  # rescale for simulations
  time <- time/cores            # rescale by number of cores
  time <- time*tincr            # rescale for number of obs
  time <- time*T1               # rescale by post intervention periods
  time <- time*I                # rescale by number of treated units
  
  time <- time/60               # Convert in minutes
  time <- ceiling(time)         # Take smallest larger integer
  time <- time*2
  
  if (name == "lasso") {
    time <- time*8
  }
  
  if (time < 60) {
    if (time < 1) {
      toprint <- "Maximum expected execution time: less than a minute.\n"
    } else if (time == 1) {
      toprint <- paste0("Maximum expected execution time: ",time," minute.\n")
    } else {
      toprint <- paste0("Maximum expected execution time: ",time," minutes.\n")
    }
  } else {
    hours <- floor(time/60)
    if (hours == 1) {
      toprint <- paste0("Maximum expected execution time: more than ",hours," hour.\n")
    } else {
      toprint <- paste0("Maximum expected execution time: more than ",hours," hours.\n")
    }
  }
  
  cat(toprint)
  cat("\n")
}


mat2list <- function(mat, cols = TRUE){
  # select rows
  names <- strsplit(rownames(mat), "\\.")
  rnames <- unlist(lapply(names, "[[", 1))
  tr.units <- unique(rnames)
  
  # select columns
  matlist <- list()
  if (cols == TRUE) {
    if (ncol(mat) > 1) {
      names <- strsplit(colnames(mat), "\\.")
      cnames <- unlist(lapply(names, "[[", 1))
      for (tr in tr.units) {
        matlist[[tr]] <- mat[rnames == tr, cnames == tr, drop=FALSE]
      }
    } else if (ncol(mat) == 1) {
      for (tr in tr.units) {
        matlist[[tr]] <- mat[rnames == tr, 1, drop=FALSE]
      }
    } else {
      for (tr in tr.units) {
        matlist[[tr]] <- mat[rnames == tr, 0, drop=FALSE]
      }
    }
  } else if (cols == FALSE) {
    for (tr in tr.units) {
      matlist[[tr]] <- mat[rnames == tr, , drop=FALSE]
    }
  }
  
  return(matlist)
}

ci2df <- function(CI, type) {
  names <- strsplit(rownames(CI), "\\.")
  df <- data.frame(ID = unlist(lapply(names, "[[", 1)),
                   Time = unlist(lapply(names, "[[", 2)),
                   lb = CI[,1], ub = CI[,2])
  
  names(df) <- c("ID","Time", paste0("lb.",type), paste0("ub.",type))
  
  return(df)
}

detectConstant <- function(x) {
  n <- nrow(x)
  col.keep <- apply(x, 2, function(j) sum(j == 1)) != n # remove double constant
  col.keep2 <- colSums(x) != 0 # remove constant other features
  x <- cbind(x[, (col.keep & col.keep2), drop = FALSE], 1)
  return(x)
}

trendRemove <- function(mat) {
  sel <- c()
  for (l in stringr::str_split(colnames(mat), "\\.")) {
    if (length(l) < 3) {
      sel <- c(sel, TRUE)
    } else {
      if (l[[3]] == "trend") {
        sel <- c(sel, FALSE)
      } else {
        sel <- c(sel, TRUE)
      }
    }
  }
  
  return(list(mat=mat[,sel,drop=FALSE], sel=sel))
}
