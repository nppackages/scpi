blockdiag <- function(I, Jtot, J, KMI, ns, slack = FALSE) {
  mat <- matrix(0, nrow = I, ncol = Jtot + KMI + ns)
  j.lb <- 1
  j.ub <- J[[1]]

  if (slack == TRUE) {
    j.lb <- j.lb + Jtot + KMI
    j.ub <- j.ub + Jtot + KMI
  }

  for (i in seq_len(I)) {
    if (i > 1) {
      j.lb <- j.ub + 1
      j.ub <- j.lb + J[[i]] - 1
    }

    mat[i, j.lb:j.ub] <- 1
  }

  return(mat)
}

blockdiagRidge <- function(Jtot, J, KMI, I) {

  mat <- matrix(0, Jtot + 2 * I, Jtot + KMI + I + 1)

  i.lb <- 1 + 2
  i.ub <- J[[1]] + 2
  j.lb <- 1
  j.ub <- J[[1]]

  for (i in seq_len(I)) {
    if (i > 1) {
      j.lb <- j.ub + 1
      j.ub <- j.lb + J[[i]] - 1
      i.lb <- i.ub + 1 + 2
      i.ub <- i.lb + J[[i]] - 1
    }

    mat[(i.lb - 2):(i.lb - 1), Jtot + KMI + i] <- c(-1, 1)
    mat[i.lb:i.ub, j.lb:j.ub] <- -diag(2, i.ub - i.lb + 1, j.ub - j.lb + 1)
  }

  return(mat)
}

ECOS_get_n_slacks <- function(w.constr, Jtot, I) {

  n_slacks <- 1

  # in lasso we add one slack per component of w to handle the abs value
  if (w.constr[["p"]] == "L1" && w.constr[["dir"]] == "<=") { # lasso
    n_slacks <- Jtot + n_slacks
  }

  # in ridge we have two hyperbolic constraints (norm and loss function)
  if (w.constr[["p"]] == "L2" && w.constr[["dir"]] == "<=") { # ridge
    n_slacks <- I + n_slacks
  }

  # L1-L2 combines ridge and simplex slack variables
  if (w.constr[["p"]] == "L1-L2") { # L1-L2
    n_slacks <- I + n_slacks
  }

  return(n_slacks)
}

ECOS_get_dims <- function(Jtot, J, KMI, w.constr, I, red) {

  if (w.constr[["p"]] == "L1" && w.constr[["dir"]] == "==") { # simplex
    dims <- list("l" = Jtot + 1, "q" = list(Jtot + KMI + 2 - red), "e" = 0)

  } else if (w.constr[["p"]] == "L1" && w.constr[["dir"]] == "<=") { # lasso
    dims <- list("l" = 1 + 2 * Jtot + I, "q" = list(Jtot + KMI + 2 - red), "e" = 0)

  } else if (w.constr[["p"]] == "L2" && w.constr[["dir"]] == "<=") { # ridge
    dims <- list("l" = 1 + I, "q" = append(lapply(J, function(i) i + 2), Jtot + KMI + 2 - red), "e" = 0)

  } else if (w.constr[["p"]] == "L1-L2") { # L1-L2
    dims <- list("l" = 1 + I + Jtot, "q" = append(lapply(J, function(i) i + 2), Jtot + KMI + 2 - red), "e" = 0)

  } else if (w.constr[["p"]] == "no norm") { # ols
    dims <- list("l" = 1, "q" = list(Jtot + KMI + 2 - red), "e" = 0)

  }

  return(dims)
}

ECOS_get_c <- function(xt, ns) {
  C <- c(xt, rep(0, ns))
  return(C)
}

ECOS_get_A <- function(J, Jtot, KMI, I, w.constr, ns) {

  if ((w.constr[["p"]] == "L1" && w.constr[["dir"]] == "==") || w.constr[["p"]] == "L1-L2") { # simplex, L1-L2

    A <- blockdiag(I, Jtot, J, KMI, ns)

  } else  { # ols, lasso, ridge

    A <- matrix(NA, 0, 0)

  }

  return(methods::as(A, "sparseMatrix"))
}

ECOS_get_b <- function(Q1, Q2, w.constr) {

  if ((w.constr[["p"]] == "L1" && w.constr[["dir"]] == "==") || w.constr[["p"]] == "L1-L2") { # simplex, L1-L2
    b <- Q1

  } else { # ols, lasso, ridge
    b <- NULL
  }

  return(b)
}

ECOS_get_G <- function(Jtot, KMI, J, I, a, Q, w.constr, ns, red, scale) {

  if (w.constr[["p"]] == "L1" && w.constr[["dir"]] == "==") { # simplex

    G <- rbind(c(a, scale),                                                              # linear part of QF
               cbind(-diag(1,Jtot), matrix(0, Jtot, KMI), matrix(0, Jtot, ns)),          # lower bounds on w
               c(rep(0, Jtot+KMI), -1),                                                  # SOC definition (||sqrt(Q)beta|| <= t)
               c(rep(0, Jtot+KMI), 1),
               cbind(-2*Q, rep(0, Jtot+KMI-red)))

  } else if (w.constr[["p"]] == "L1" && w.constr[["dir"]] == "<=") { # lasso, x = (beta, z, t)

    G <- rbind(c(a, rep(0, ns - 1), scale),                                              # linear part of QF
               cbind(diag(1, Jtot), matrix(0, Jtot, KMI), diag(1, Jtot), rep(0, Jtot)),  # z \geq -w
               cbind(-diag(1, Jtot), matrix(0, Jtot, KMI), diag(1, Jtot), rep(0, Jtot)), # z \geq w
               -blockdiag(I, Jtot, J, KMI, ns, TRUE),                                    # norm-inequality constraint
               c(rep(0, Jtot+KMI+Jtot), -1),                                             # SOC definition (||sqrt(Q)beta|| <= t)
               c(rep(0, Jtot+KMI+Jtot), 1),
               cbind(-2 * Q, matrix(0, Jtot + KMI - red, Jtot + 1)))

  } else if (w.constr[["p"]] == "L2" && w.constr[["dir"]] == "<=") { # ridge, x = (beta, s, t)

    G <- rbind(c(a, rep(0, I), scale),                                                   # linear part of QF
               cbind(matrix(0, I, Jtot + KMI), diag(1, I, I), rep(0, I)),                # s \leq Q1^2
               blockdiagRidge(Jtot, J, KMI, I),                                          # SOC definition (||w|| <= s)
               c(rep(0, Jtot + KMI), rep(0, I), -1),                                     # SOC definition (||sqrt(Q)beta|| <= t)
               c(rep(0, Jtot + KMI), rep(0, I), 1),
               cbind(-2 * Q, matrix(0, Jtot + KMI - red, I + 1)))

  } else if (w.constr[["p"]] == "L1-L2") { # L1-L2, x = (beta, s, t)

    G <- rbind(c(a, rep(0, ns-1), scale),                                                # linear part of QF
               cbind(-diag(1,Jtot), matrix(0, Jtot, KMI), matrix(0, Jtot, ns)),          # lower bounds on w
               cbind(matrix(0, I, Jtot+KMI), diag(1, I, I), rep(0, I)),                  # s \leq Q2^2
               blockdiagRidge(Jtot, J, KMI, I),                                          # SOC definition (||w||_2 <= s)
               c(rep(0, Jtot+KMI), rep(0, I), -1),                                       # SOC definition (||sqrt(Q)beta||_2 <= t)
               c(rep(0, Jtot+KMI), rep(0, I), 1),
               cbind(-2*Q, matrix(0, Jtot + KMI - red, I + 1)))

  } else if (w.constr[["p"]] == "no norm") { # ols

    G <- rbind(c(a, scale),                                                              # linear part of QF
               c(rep(0, Jtot+KMI), -1),                                                  # SOC definition (||sqrt(Q)beta|| <= t)
               c(rep(0, Jtot+KMI), 1),
               cbind(-2 * Q, rep(0, Jtot + KMI - red)))

  }

  return(methods::as(G, "sparseMatrix"))
}


ECOS_get_h <- function(d, lb, J, Jtot, KMI, I, w.constr, Q1, Q2, red) {

  if (w.constr[["p"]] == "L1" && w.constr[["dir"]] == "==") { # simplex

    h <- c(-d,                         # linear part of QF
           -lb,                        # lower bounds of w
           1, 1, rep(0,Jtot+KMI-red))  # SOC definition

  } else if (w.constr[["p"]] == "L1" && w.constr[["dir"]] == "<=") { # lasso

    h <- c(-d,                         # linear part of QF 
           rep(0, 2*Jtot),             # abs(w) <= z
           Q1,                         # norm-inequality constraints
           1, 1, rep(0,Jtot+KMI-red))  # SOC definition

  } else if (w.constr[["p"]] == "L2" && w.constr[["dir"]] == "<=") { # ridge

    aux <- unlist(lapply(1:I, function(x) c(1,1,rep(0,J[[x]]))))

    h <- c(-d,                         # linear part of QF
           Q1^2,                       # s <= Q1^2
           aux,                        # SOC definition (||w|| <= s)
           1, 1, rep(0,Jtot+KMI-red))  # SOC definition (||sqrt(Q)beta|| <= t)

  } else if (w.constr[["p"]] == "L1-L2") { # L1-L2

    aux <- unlist(lapply(1:I, function(x) c(1, 1, rep(0, J[[x]]))))

    h <- c(-d,                         # linear part of QF
           -lb,                        # lower bounds of w
           Q2^2,                       # s <= Q2^2
           aux,                        # SOC definition (||w||_2 <= s)
           1, 1, rep(0,Jtot+KMI-red))  # SOC definition (||sqrt(Q)beta|| <= t)

  } else if (w.constr[["p"]] == "no norm") { # ols

    h <- c(-d,                         # linear part of QF
           1, 1, rep(0,Jtot+KMI-red))  # SOC definition

  }

  return(h)
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
      w.constr[["Q"]]      <- max(min(Qfeat, na.rm = TRUE), .5)
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
        Af <- A[feature.id == feat, , drop = FALSE]
        Zf <- Z[feature.id == feat, , drop = FALSE]
        Vf <- V[feature.id == feat, feature.id == feat, drop = FALSE]

        if (nrow(Af) >= 5) {
          QQ <- tryCatch({
            aux <- shrinkage.EST("ridge", Af, Zf, Vf, J, KM)
            Q2 <- aux$Q
          }, warning = {}, error = {}, finally = {})
          Qfeat <- c(Qfeat, QQ)
        }
      }

      if (is.null(Qfeat)) Qfeat <- shrinkage.EST("ridge", A, Z, V, J, KM)$Q
      w.constr[["Q2"]]      <-  max(min(Qfeat, na.rm = TRUE), .5)
      w.constr[["lambda"]] <- aux$lambda
    }

    w.constr <- list(lb     = 0,
                     p      = "L1-L2",
                     dir    = "==/<=",
                     Q      = 1,
                     Q2     = w.constr[["Q2"]],
                     name   = "L1-L2",
                     lambda = w.constr[["lambda"]])

  } else {

    # if constraint is entirely user specified just check everything is fine
    if (!(all(c('p', 'dir', 'Q', 'lb') %in% names(w.constr)))) {
      stop("If 'name' is not specified, w.constr should be a list whose elements 
            must be named 'p','dir','Q','lb'.")
    }

    if (!(w.constr[["p"]] %in% c("no norm", "L1", "L2", "L1-L2"))) {
      stop("Specify either p = 'no norm' (no constraint on the norm of weights),
          p = 'L1' (L1-norm), p = 'L2' (L2-norm)")
    } else if (w.constr[["p"]] == "no norm") {
      w.constr[["dir"]] <- "NULL"
    }

    if (!(w.constr[["dir"]] %in% c("<=", "==", "==/<=", "NULL"))) {
      stop("Specify either dir = '<=' (inequality constraint on the norm of the weights)
            or dir = '==' (equality constraint on the norm of the weights) or dir = 'NULL'
            in case you don't want to specify a constraint on the norm of the weights.")
    }

    if (!(w.constr[["lb"]] == 0 || w.constr[["lb"]] == -Inf)) {
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
    lambd   <- sig.wls^2 * (J + KM) / sum(wls$coef^2, na.rm = TRUE)           # rule of thumb for lambda (Hoerl et al, 1975)
    Q       <- sqrt(sum(wls$coef^2, na.rm = TRUE)) / (1 + lambd)              # convert lambda into Q

    if (is.nan(Q) || (nrow(Z) <= ncol(Z) + 10)) { # reduce dimensionality of the problem if more params than obs
      lasso.cols <- b.est(A, Z, J, KM, list(dir = "<=", lb = -Inf, p = "L1", Q = 1), V)
      active.cols <- abs(lasso.cols) > 1e-8
      if (sum(active.cols) >= (max(nrow(A) - 10, 2)) ) {
        active.cols <-  rank(-abs(lasso.cols)) <= max(nrow(A) - 10, 2)
      }
      Z.sel <- Z[, active.cols, drop = FALSE]
      wls     <- lm(A ~ Z.sel - 1, weights = diag(V))
      sig.wls <- sigma(wls)
      lambd   <- sig.wls^2 * (ncol(Z.sel) + KM) / sum(wls$coef^2, na.rm = TRUE)     # rule of thumb for lambda (Hoerl et al, 1975)
      Q       <- sqrt(sum(wls$coef^2, na.rm = TRUE)) / (1 + lambd)                  # convert lambda into Q
    }
  }

  return(list(Q = Q, lambda = lambd))
}

# Auxiliary function that solves the (un)constrained problem to estimate b
# depending on the desired method
b.est <- function(A, Z, J, KM, w.constr, V, CVXR.solver = "ECOS") {

  dire <- w.constr[["dir"]]
  lb   <- w.constr[["lb"]]
  p    <- w.constr[["p"]]
  QQ   <- w.constr[["Q"]]

  if (p == "L1-L2") {
    Q2 <- w.constr[["Q2"]]
  } else {
    Q2 <- NULL
  }

  x <- CVXR::Variable(J + KM)
  objective <- CVXR::Minimize(CVXR::quad_form(A - Z %*% x, V))

  if (p == "no norm") { # least squares
    constraints <- list()

  } else if (p == "L1") {
    if (dire == "==") { # simplex
      constraints <- list(CVXR::sum_entries(x[1:J]) == QQ, x[1:J] >= lb)
    } else if (dire == "<=") { # lasso
      constraints <- list(CVXR::norm1(x[1:J]) <= QQ, x[1:J] >= lb)
    }

  } else if (p == "L2") { # ridge
    if (dire == "==") {
      constraints <- list(CVXR::sum_squares(x[1:J]) == CVXR::power(QQ, 2))
    } else if (dire == "<=") {
      constraints <- list(CVXR::sum_squares(x[1:J]) <= CVXR::power(QQ, 2))
    }

  } else if (p == "L1-L2") {
    constraints <- list(CVXR::sum_entries(x[1:J]) == QQ,
                        CVXR::power(CVXR::cvxr_norm(x[1:J], 2), 2) <= CVXR::power(Q2, 2),
                        x[1:J] >= lb)
  }

  prob        <- CVXR::Problem(objective, constraints)
  sol         <- CVXR::solve(prob, solver=CVXR.solver)
  
  b <- sol$getValue(x)
  alert <- !(sol$status %in% c("optimal", "optimal_inaccurate"))

  if (alert == TRUE) {
    stop(paste0("Estimation algorithm not converged! The algorithm returned the value:",
                sol$status, ". To check to what errors it corresponds go to 
               'https://cvxr.rbind.io/cvxr_examples/cvxr_gentle-intro/'. Typically, this occurs
                because the problem is badly-scaled. If so, scaling the data fixes the issue. Another
                fix could be changing the algorithm via the option 'solver'. Check your available options
                using CVXR::installed_solvers() and consult 
                'https://cvxr.rbind.io/cvxr_examples/cvxr_using-other-solvers/'"),
         call. = FALSE)
  }

  b <- b[, 1, drop = TRUE]
  names(b) <- colnames(Z)

  return(b)
}

# Auxiliary function that solves the (un)constrained problem to estimate b
# depending on the desired method - Multiple treated units case
b.est.multi <- function(A, Z, J, KMI, I, w.constr, V, CVXR.solver="ECOS") {

  # The constraint is symmetric in the shape across treated units (J, KM, Q might change)
  dire  <- w.constr[[1]]$dir
  lb    <- w.constr[[1]]$lb
  p     <- w.constr[[1]]$p
  QQ    <- unlist(lapply(w.constr, function(x) x$Q))

  if (p == "L1-L2") {
    Q2 <- unlist(lapply(w.constr, function(x) x$Q2))
  }

  Jtot <- sum(unlist(J))

  x <- CVXR::Variable(Jtot + KMI)

  objective   <- CVXR::Minimize(CVXR::quad_form(A - Z %*% x, V))

  if (lb != -Inf) {
    constraints <- list(x[1:Jtot] >= lb)
  } else {
    constraints <- list()
  }

  j.lb <- 1
  for (i in seq_len(I)) {
    j.ub <- j.lb + J[[i]] - 1

    if (p == "L1") {
      if (dire == "==") { # simplex
        constraints <- append(constraints, list(CVXR::sum_entries(x[j.lb:j.ub]) == QQ[i]))
      } else if (dire == "<=") { # lasso
        constraints <- append(constraints, list(CVXR::norm1(x[j.lb:j.ub]) <= QQ[i]))
      }

    } else if (p == "L2") { # ridge
        constraints <- append(constraints, list(CVXR::sum_squares(x[j.lb:j.ub]) <= QQ[i]^2))

    } else if (p == "L1-L2") {
      constraints <- append(constraints, list(CVXR::sum_entries(x[j.lb:j.ub]) == QQ[i],
                          CVXR::power(CVXR::cvxr_norm(x[j.lb:j.ub], 2), 2) <= CVXR::power(Q2[i], 2)))
    }

    j.lb <- j.ub + 1
  }

  prob        <- CVXR::Problem(objective, constraints)
  sol         <- CVXR::solve(prob, solver=CVXR.solver)

  alert <- !(sol$status %in% c("optimal", "optimal_inaccurate"))
  
  if (alert == TRUE) {
    stop(paste0("Estimation algorithm not converged! The algorithm returned the value:",
                sol$status, ". To check to what errors it corresponds go to 
               'https://cvxr.rbind.io/cvxr_examples/cvxr_gentle-intro/'. Typically, this occurs
                because the problem is badly-scaled. If so, scaling the data fixes the issue. Another
                fix could be changing the algorithm via the option 'solver'. Check your available options
                using CVXR::installed_solvers() and consult 
                'https://cvxr.rbind.io/cvxr_examples/cvxr_using-other-solvers/'"),
         call. = FALSE)
  }

  rownames(b)  <- colnames(Z)

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
    VV <- kronecker(ones %*% t(ones), eye) # structure if T0*M was balanced across treated unit

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
        BB       <- B[feature.id == feature, ]
        B.diff   <- rbind(B.diff, BB - dplyr::lag(BB))
      }
      u.des.0 <- cbind(B.diff, C)[, index, drop = FALSE]    # combine with C

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
        BB       <- B[feature.id == feature, ]
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
    colnames(B.lag) <- rep(paste0(name.tr, ".lag"), ncol(B.lag))
    u.des.0 <- cbind(u.des.0, B.lag[, index.w, drop = FALSE])
  }

  # If user provided check compatibility of the matrix and overwrite what has been created
  if (is.null(u.design) == FALSE) {
    if (is.matrix(u.design) == FALSE) {
      stop("The object u.design should be a matrix!!")
    }

    if (nrow(u.design) != nrow(res)) {
      stop(paste("The matrix u.design has", nrow(u.design), "rows when", nrow(res),
                 "where expected!"))
    }
    u.des.0 <- u.design
  }

  return(list(u.des.0 = u.des.0))
}


e.des.prep <- function(B, C, P, e.order, e.lags, res, sc.pred, Y.donors, out.feat, features,
                       J, index, index.w, coig.data, T0, T1, constant, e.design, P.diff.pre,
                       effect, I, class.type) {

  aux <- trendRemove(P)
  C <- trendRemove(C)$mat
  index <- index[aux$sel]
  P <- aux$mat

  if (!is.null(P.diff.pre)) P.diff.pre <- trendRemove(as.matrix(P.diff.pre))$mat

  if (out.feat == FALSE) {
    # If the outcome variable is not among the features we need to create a
    # proper vector of residuals. Further, we force the predictors to be
    # the outcome variable of the donors
    if (class.type == "scpi_data") { # just one treated unit
      e.res <- sc.pred$data$Y.pre - sc.pred$est.results$Y.pre.fit
    } else { # need to extract data from the specific treated unit
      tr.unit <- stringr::str_split(rownames(res[1,,drop=FALSE]), "\\.")[[1]][[1]]
      sel <- sapply(stringr::str_split(rownames(sc.pred$data$Y.pre), "\\."), "[[", 1) == tr.unit
      e.res <- sc.pred$data$Y.pre[sel, 1, drop=FALSE] - sc.pred$est.results$Y.pre.fit[sel, 1, drop=FALSE]
    }

    if (coig.data == TRUE) {
      e.des.0 <- apply(Y.donors, 2, function(x) x - dplyr::lag(x))[, index.w, drop=FALSE]

      if (effect == "time") {
        P.first <- (P[1, ]*I - Y.donors[T0[1], ])/I
      } else {
        P.first <- P[1, ] - Y.donors[T0[1], ]
      }
      P.diff  <- rbind(P.first, apply(P, 2, diff))[, index, drop = FALSE]
      e.des.1 <- P.diff

    } else {
      e.des.0  <- Y.donors[, index.w, drop=FALSE]
      e.des.1  <- P[, index, drop = FALSE]
    }
    
  } else if (out.feat == TRUE) {    # outcome variable is among features
    e.res <- res[1:T0[1], , drop = FALSE]

    ## Construct the polynomial terms in B (e.des.0) and P (e.des.1)
    if (e.order == 0) {

      e.des.0 <- as.matrix(rep(1, T0[1]))
      if (effect == "time") {
        e.des.1 <- as.matrix(rep(1/I, T1))
      } else {
        e.des.1 <- as.matrix(rep(1, T1))
      }

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
        if (effect == "time") {
          P.first <- c((P[1, (1:J), drop = FALSE]*I - B[feature.id == features[1], , drop = FALSE][T0[1], ]),
                       P[1, -(1:J), drop = FALSE]*I)/I

        } else {
          P.first <- c((P[1, (1:J), drop = FALSE] - B[feature.id == features[1], , drop = FALSE][T0[1], ]),
                       P[1, -(1:J), drop = FALSE])
        }

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
        if (effect == "time") {
          e.des.1 <- cbind(e.des.1, rep(1/I, nrow(e.des.1)))
        } else {
          e.des.1 <- cbind(e.des.1, rep(1, nrow(e.des.1)))
        }
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


insampleUncertaintyGet <- function(Z.na, V.na, P.na, beta, Sigma.root, J, KMI, I,
                                   w.constr.inf, Q.star, Q2.star, lb, TT, sims, cores, verbose, 
                                   w.lb.est, w.ub.est) {

  Q <- t(Z.na) %*% V.na %*% Z.na / TT
  colnames(Q) <- colnames(Z.na)

  jj <- nrow(P.na)

  # optimizing options

  Jtot <- sum(unlist(J))

  iters <- round(sims / 10)
  perc  <- 0

  # simulate
  ns <- ECOS_get_n_slacks(w.constr.inf, Jtot, I)
  decomp <- matRegularize(Q)
  Qreg <- decomp$Qreg
  scale <- decomp$scale
  dimred <- nrow(Q) - nrow(Qreg)

  data <- list()
  data[["dims"]] <- ECOS_get_dims(Jtot, J, KMI, w.constr.inf, I, dimred)
  data[["A"]] <- ECOS_get_A(J, Jtot, KMI, I, w.constr.inf, ns)
  data[["b"]] <- ECOS_get_b(Q.star, Q2.star, w.constr.inf)

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
      G <- Sigma.root %*% zeta

      a <- -2 * G - 2 * Q %*% beta
      d <- 2 * sum(G * beta) + sum(beta * (Q %*% beta))
      if (methods::is(Q, "dgeMatrix") == TRUE) a <- as.matrix(a)
 
      data[["G"]] <- ECOS_get_G(Jtot, KMI, J, I, a, Qreg, w.constr.inf, ns, dimred, scale)
      data[["h"]] <- ECOS_get_h(d, lb, J, Jtot, KMI, I, w.constr.inf, Q.star, Q2.star, dimred)

      for (hor in seq_len(jj)) {
        xt <- P.na[hor, ]
        # minimization
        if (w.lb.est == TRUE) {
          data[["c"]] <- ECOS_get_c(-xt, ns)

          solver_output <- ECOSolveR::ECOS_csolve(c = data[["c"]],
                                                  G = data[["G"]],
                                                  h = data[["h"]],
                                                  dims = data[["dims"]],
                                                  A = data[["A"]],
                                                  b = data[["b"]])

          if (!(solver_output$infostring %in% c("Optimal solution found", "Close to optimal solution found"))) {
            lb.f <- NA
          } else {
            xx <- solver_output$x[1:(Jtot + KMI)]
            lb.f <- -sum(xt * (xx - beta))
          }
        } else {
          lb.f <- NA
        }

        # maximization
        if (w.ub.est == TRUE) {
          data[["c"]] <- ECOS_get_c(xt, ns)

          solver_output <- ECOSolveR::ECOS_csolve(c = data[["c"]],
                                                  G = data[["G"]],
                                                  h = data[["h"]],
                                                  dims = data[["dims"]],
                                                  A = data[["A"]],
                                                  b = data[["b"]])

          if (!(solver_output$infostring %in% c("Optimal solution found", "Close to optimal solution found"))) {
            ub.f <- NA
          } else {
            xx <- solver_output$x[1:(Jtot+KMI)]
            ub.f <- -sum(xt*(xx - beta))
          }     
        } else {
          ub.f <- NA
        }

        vsig[sim, hor] <- lb.f
        vsig[sim, hor + nrow(P.na)] <- ub.f
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
                             .packages = c('ECOSolveR', 'Matrix', 'methods'),
                             .export   = c('ECOS_get_G', 'ECOS_get_h', 'ECOS_get_h', 'ECOS_get_c',
                                           'blockdiag', 'blockdiagRidge', 'sqrtm'),
                             .combine  = rbind,
                             .options.snow = opts) %dopar% {

                               ub.sim <- lb.sim <- c()

                               zeta    <- rnorm(length(beta))
                               G       <- Sigma.root %*% zeta
                               a <- -2 * G - 2 * Q %*% beta
                               d <- 2 * sum(G * beta) + sum(beta * (Q %*% beta))
                               if (methods::is(Q, "dgeMatrix") == TRUE) a <- as.matrix(a)
                               
                               data[["G"]] <- ECOS_get_G(Jtot, KMI, J, I, a, Qreg, w.constr.inf, ns, dimred, scale)
                               data[["h"]] <- ECOS_get_h(d, lb, J, Jtot, KMI, I, w.constr.inf, Q.star, Q2.star, dimred)

                               for (hor in seq_len(jj)) {
                                 xt <- P.na[hor, ]

                                 # minimization
                                 if (w.lb.est == TRUE) {
                                   data[["c"]] <- ECOS_get_c(-xt, ns)

                                   solver_output <- ECOSolveR::ECOS_csolve(c = data[["c"]],
                                                                           G = data[["G"]],
                                                                           h = data[["h"]],
                                                                           dims = data[["dims"]],
                                                                           A = data[["A"]],
                                                                           b = data[["b"]])

                                   if (!(solver_output$infostring %in% c("Optimal solution found", "Close to optimal solution found"))) {
                                     lb.f <- NA
                                   } else {
                                     xx <- solver_output$x[1:(Jtot + KMI)]
                                     lb.f <- -sum(xt * (xx - beta))
                                   }
                                 } else {
                                   lb.f <- NA
                                 }

                                 # maximization
                                 if (w.ub.est == TRUE) {
                                   data[["c"]] <- ECOS_get_c(xt, ns)

                                   solver_output <- ECOSolveR::ECOS_csolve(c = data[["c"]],
                                                                           G = data[["G"]],
                                                                           h = data[["h"]],
                                                                           dims = data[["dims"]],
                                                                           A = data[["A"]],
                                                                           b = data[["b"]])

                                   if (!(solver_output$infostring %in% c("Optimal solution found", "Close to optimal solution found"))) {
                                     ub.f <- NA
                                   } else {
                                     xx <- solver_output$x[1:(Jtot+KMI)]
                                     ub.f <- -sum(xt * (xx - beta))
                                   }
                                 } else {
                                   ub.f <- NA
                                 }

                                 lb.sim      <- append(lb.sim, lb.f)
                                 ub.sim      <- append(ub.sim, ub.f)

                               }

                               c(lb.sim, ub.sim)
                             }

    parallel::stopCluster(cl)
  }

  return(vsig)
}

# Prediction interval, for e
scpi.out <- function(res, x, eval, e.method, alpha, e.lb.est, e.ub.est, verbose, effect, out.feat) {

  neval <- nrow(eval)
  e.1 <- e.2 <- lb <- ub <- NA

  if (e.lb.est == TRUE || e.ub.est == TRUE) {
    if (e.method == "gaussian") {
      x.more   <- rbind(eval, x)
      fit      <- predict(y = res, x = x, eval = x.more, type = "lm")
      e.mean   <- fit[1:neval]
      res.fit  <- fit[-(1:neval)]

      if (effect == "time") {
        labels <- unlist(purrr::map(stringr::str_split(rownames(fit)[1:neval], "\\."), 2))
        e.mean <- aggregateUnits(xx = e.mean, labels = labels)
      }

      var.pred <- predict(y = log((res - res.fit)^2), x = x, eval = x.more, type = "lm")
      if (effect == "time") {
        var.pred <- aggregateUnits(xx = var.pred[1:neval], labels = labels)
      } else {
        var.pred <- var.pred[1:neval]
      }
      e.sig2   <- exp(var.pred)

      q.pred <- predict(y = res - res.fit, x = x, eval = x.more, type = "qreg", tau = c(0.25, 0.75))
      if (effect == "time") {
        q3.pred <- aggregateUnits(xx = q.pred[1:neval, 2], labels = labels)
        q1.pred <- aggregateUnits(xx = q.pred[1:neval, 1], labels = labels)
      } else {
        q3.pred <- q.pred[1:neval, 2]
        q1.pred <- q.pred[1:neval, 1]
      }
      IQ.pred <- q3.pred - q1.pred
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

      if (effect == "time") {
        labels <- unlist(purrr::map(stringr::str_split(rownames(fit)[1:neval], "\\."), 2))
        e.mean <- aggregateUnits(xx=e.mean, labels=labels)
      }

      var.pred <- predict(y=log((res-res.fit)^2), x=x, eval=x.more, type="lm")
      res.var <- var.pred[-(1:neval)]
      if (effect == "time") {
        var.pred <- aggregateUnits(xx=var.pred[1:neval], labels=labels)
      } else {
        var.pred <- var.pred[1:neval]
      }
      e.sig    <- sqrt(exp(var.pred))
      res.st   <- (res-res.fit)/sqrt(exp(res.var))

      q.pred <- predict(y=res-res.fit, x=x, eval=x.more, type="qreg", tau = c(0.25,0.75))
      if (effect == "time") {
        q3.pred <- aggregateUnits(xx=q.pred[1:neval,2], labels=labels)
        q1.pred <- aggregateUnits(xx=q.pred[1:neval,1], labels=labels)
      } else {
        q3.pred <- q.pred[1:neval,2]
        q1.pred <- q.pred[1:neval,1]
      }
      IQ.pred <- q3.pred - q1.pred
      IQ.pred <- abs(IQ.pred)
      e.sig <- apply(cbind(e.sig, IQ.pred/1.34), 1, min)

      lb <- e.mean + e.sig * quantile(res.st, alpha)
      ub <- e.mean + e.sig * quantile(res.st, 1-alpha)

      # save mean and variance of u, only for sensitivity analysis
      e.1 <- e.mean
      e.2 <- e.sig^2

    } else if (e.method == "qreg") {
      e.pred  <- predict(y=res, x=x, eval=eval, type="qreg", tau=c(alpha, 1-alpha), verbose = verbose)

      if (effect == "time") {
        labels <- unlist(purrr::map(stringr::str_split(rownames(e.pred)[1:neval], "\\."), 2))
        e.pred.lb <- aggregateUnits(e.pred[,1], labels)
        e.pred.ub <- aggregateUnits(e.pred[,2], labels)
      } else {
        e.pred.lb <- e.pred[,1]
        e.pred.ub <- e.pred[,2]
      }

      lb <- e.pred.lb
      ub <- e.pred.ub

    }
  }

  # if model is heavily misspecified just give up on out-of-sample uncertainty
  if (out.feat[[1]] == FALSE) {
    if (any(lb > 0)) lb <- rep(min(lb), length(lb))
    if (any(ub < 0)) ub <- rep(max(ub), length(ub))
  }
  
  return(list(lb = lb, ub = ub, e.1 = e.1, e.2 = e.2))
}


simultaneousPredGet <- function(vsig, T1, T1.tot, I, u.alpha, e.alpha, e.res.na, e.des.0.na, e.des.1,
                                w.lb.est, w.ub.est, w.bounds, w.name, effect, out.feat) {

  vsigUB <- vsig[, (T1.tot + 1):(2 * T1.tot), drop = FALSE]
  vsigLB <- vsig[, 1:T1.tot, drop = FALSE]
  
  pi.e   <- scpi.out(res = e.res.na, x = e.des.0.na, eval = e.des.1, 
                     e.method = "gaussian", alpha = e.alpha/2, e.lb.est = TRUE,
                     e.ub.est =  TRUE, effect = effect, out.feat = out.feat)

  w.lb.joint <- w.ub.joint <- c()
  
  j.min <- 1
  
  for (i in seq_len(I)) {
    j.max <- T1[[i]] + j.min - 1

    lb.joint <- quantile(apply(vsigLB[, j.min:j.max, drop = FALSE], 1, min, na.rm = TRUE), probs = u.alpha/2)
    ub.joint <- quantile(apply(vsigUB[, j.min:j.max, drop = FALSE], 1, max, na.rm = TRUE), probs = (1-u.alpha/2))

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

epskappaGet <- function(P, rho.vec, beta, I, effect, joint = FALSE) {
  beta.list <- mat2list(beta)
  
  if (effect == "time") {
    rho.vec <- mean(rho.vec)
    I <- 1
    P.list <- list(P)
    beta.list <- list(beta)
  } else {
    P.list <- mat2list(P)
  }
  
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
  
  if (u.sigma == "HC0") { # White (1980)
    vc <- 1
  }
  
  else if (u.sigma == "HC1") { # MacKinnon and White (1985)
    vc <- TT/(TT-df)
  }
  
  else if (u.sigma == "HC2") { # MacKinnon and White (1985)
    PP <- Z %*% Matrix::solve(t(Z)%*% V %*% Z) %*% t(Z) %*% V
    vc <- 1/(1-diag(PP))
  }
  
  else if (u.sigma == "HC3") { # Davidson and MacKinnon (1993)
    PP <- Z %*% Matrix::solve(t(Z)%*% V %*% Z) %*% t(Z) %*% V
    vc <- 1/(1-diag(PP))^2
  }
  
  else if (u.sigma == "HC4") { # Cribari-Neto (2004)
    PP <- Z %*% Matrix::solve(t(Z)%*% V %*% Z) %*% t(Z) %*% V
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
    denomCheck(sigma.bj)
    CC       <- sigma.u/sigma.bj

  } else if (rho == "type-2"){
    sigma.u   <- sqrt(mean((res-mean(res))^2))
    sigma.bj2 <- min(apply(B, 2, var))
    denomCheck(sigma.bj2)
    sigma.bj  <- max(apply(B, 2, sd))
    CC        <- sigma.bj*sigma.u/sigma.bj2

  } else if (rho == "type-3"){
    sigma.bj2 <- min(apply(B, 2, var))
    denomCheck(sigma.bj2)
    sigma.bju <- max(apply(B, 2, function(bj) cov(bj, res)))
    CC        <- sigma.bju/sigma.bj2
  }

  if (coig.data == TRUE) { # cointegration
    c <- 1
  } else {        # iid or ar
    c <- 0.5
  }

  rho <- (CC * (log(T0.tot))^c) / (sqrt(T0.tot))

  if (is.null(rho.max) == FALSE)  rho <- min(rho, rho.max)

  return(rho)
}

denomCheck <- function(den) {
  if (den == 0) {
    stop("One of your donors has no variation in the pre-treatment period!")
  }
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

# executionTime <- function(T0, J, I, T1, sims, cores, name){
#   Ttot <- sum(T0)
#   tincr <- Ttot/1000
#   coefsJ <- c(-0.54755616,0.09985644)
#   
#   time <- sum(c(1,J)*coefsJ)    # Time for J donors, T0 = 1000, sims = 10
#   time <- ceiling(time)*sims/10  # rescale for simulations
#   time <- time/cores            # rescale by number of cores
#   time <- time*tincr            # rescale for number of obs
#   time <- time*T1               # rescale by post intervention periods
#   time <- time*I                # rescale by number of treated units
#   
#   time <- time/60               # Convert in minutes
#   time <- ceiling(time)         # Take smallest larger integer
#   time <- time*2
#   
#   if (name == "lasso") {
#     time <- time*8
#   }
#   
#   if (time < 60) {
#     if (time < 1) {
#       toprint <- "Maximum expected execution time: less than a minute.\n"
#     } else if (time == 1) {
#       toprint <- paste0("Maximum expected execution time: ",time," minute.\n")
#     } else {
#       toprint <- paste0("Maximum expected execution time: ",time," minutes.\n")
#     }
#   } else {
#     hours <- floor(time/60)
#     if (hours == 1) {
#       toprint <- paste0("Maximum expected execution time: more than ",hours," hour.\n")
#     } else {
#       toprint <- paste0("Maximum expected execution time: more than ",hours," hours.\n")
#     }
#   }
#   
#   cat(toprint)
#   cat("\n")
# }


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

detectConstant <- function(x, scale.x = 1) {
  x <- x*scale.x
  n <- nrow(na.omit(x))
  col.keep <- apply(x, 2, function(j) sum(j == 1, na.rm = TRUE)) != n # remove double constant
  col.keep2 <- colSums(x, na.rm = TRUE) != 0 # remove constant other features
  x <- cbind(x[, (col.keep & col.keep2), drop = FALSE], 1)
  x <- x/scale.x
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

aggregateUnits <- function(xx, labels) {
  xx.df <- data.frame(xx = xx, id = labels)
  x.mean <- aggregate(xx.df[,"xx"], by=list(unit=xx.df$id), FUN=mean, na.rm=TRUE)
  x.mean <- x.mean[order(as.numeric(x.mean$unit)),]$x
  return(x.mean)
}

# regularizer for symmetric positive definite matrices
matRegularize <- function(P, cond = NA) {
  eig <- eigen(P, only.values = FALSE)
  w <- eig$values
  V <- eig$vectors
  
  if(is.na(cond))
    cond <- 1e6 * .Machine$double.eps   # All real numbers are stored as double precision in R
  
  scale <- max(abs(w))
  
  if(scale < cond) {
    return(list(scale = 0, M1 = V[,FALSE], M2 = V[,FALSE]))
  }
  
  w_scaled <- w / scale
  maskp <- w_scaled > cond
  maskn <- w_scaled < -cond
    
  if(any(maskp) && any(maskn)) {
    warning("Forming a non-convex expression QuadForm(x, indefinite)")
  }
  
  if(sum(maskp) <= 1) {
    M1 <- as.matrix(V[,maskp] * sqrt(w_scaled[maskp]))
  } else {
    M1 <- V[,maskp] %*% diag(sqrt(w_scaled[maskp]))
  }
  
  if(sum(maskn) <= 1){
    M2 <- as.matrix(V[,maskn]) * sqrt(-w_scaled[maskn])
  } else {
    M2 <- V[,maskn] %*% diag(sqrt(-w_scaled[maskn]))
  }
  
  if(!is.null(dim(M1)) && prod(dim(M1)) > 0) {
    Qreg <- t(M1)
  } else if (!is.null(dim(M2)) && prod(dim(M2)) > 0) {
    scale <- -scale
    Qreg <- t(M2)
  }    
  return(list(scale = scale, Qreg = Qreg))
}


outcomeGet <- function(Y.pre.fit, Y.post.fit, Y.df, units.est, treated.units,
                       plot.type, anticipation, period.post) {

  # shortcut to avoid "no visible binding for global variable 'X' when checking the package
  Treatment <- ID <- Time <- Tdate <- Type <- tstd <- NULL
      
  synth.mat <- rbind(Y.pre.fit, Y.post.fit)
  names(Y.df) <- c(names(Y.df)[1:3], "Actual")
  res.df <- subset(Y.df, ID %in% units.est)
  
  if (plot.type == "unit") {
    Y.actual.pre <- subset(res.df, Treatment == 0)
    Y.actual.post <- subset(res.df, Treatment == 1)
    Y.actual.post.agg <- aggregate(Actual ~ ID, data = Y.actual.post, mean)
    Y.actual.post.agg$Treatment <- 1
    names <- strsplit(rownames(Y.post.fit), "\\.")
    Y.actual.post.agg$Time <- as.numeric(unlist(lapply(names, "[[", 2)))
    Y.actual <- rbind(Y.actual.pre, Y.actual.post.agg)
    
  } else {
    Y.actual <- res.df # dataframe
  }
  
  treated.periods <- subset(Y.actual, Treatment == 1, select = c(Time, ID)) # post treatment period for each treated unit
  treated.reception <- aggregate(Time ~ ID, data = treated.periods, min)
  names(treated.reception) <- c("ID", "Tdate")
  treated.reception$Tdate <- as.numeric(treated.reception$Tdate) - anticipation - 1 / 2
  treated.reception <- subset(treated.reception, ID %in% units.est)
  
  if (plot.type == "time") {
    res.df <- merge(res.df, treated.reception, by = "ID")
    Y.actual.pre <- subset(res.df, Time < Tdate)
    Y.actual.post <- subset(res.df, Time > Tdate)
    Y.actual.pre$Tdate <- Y.actual.pre$Tdate + 1/2
    Y.actual.post$Tdate <- Y.actual.post$Tdate + 1/2
    Y.actual.pre$tstd <- Y.actual.pre$Time - Y.actual.pre$Tdate
    Y.actual.post$tstd <- Y.actual.post$Time - Y.actual.post$Tdate
    
    names <- strsplit(rownames(synth.mat), "\\.")
    unit <- unlist(lapply(names, "[[", 1))
    no.agg <- unit %in% treated.units
    time <- unlist(lapply(names[1:sum(no.agg)], "[[", 2))
    synth.pre <- data.frame(ID = unit[no.agg == TRUE],
                            Synthetic = synth.mat[no.agg == TRUE, 1],
                            Time = time)
    
    Y.pre <- merge(Y.actual.pre, synth.pre, by=c("ID", "Time"))
    max.pre <- max(aggregate(tstd ~ ID, data = Y.pre, min)$tstd)
    min.post <- min(unlist(lapply(period.post, length))) - 1
    
    Y.pre.agg <- aggregate(x = Y.pre[c("Actual", "Synthetic")],
                           by = Y.pre[c("tstd")],
                           FUN = mean, na.rm = TRUE)
    names(Y.pre.agg) <- c("Time", "Actual", "Synthetic")
    Y.pre.agg <- subset(Y.pre.agg, Time >= max.pre)
    
    Y.post.agg <- aggregate(x = Y.actual.post[c("Actual")],
                            by = Y.actual.post[c("tstd")],
                            FUN = mean, na.rm = TRUE)
    
    Y.post.agg <- subset(Y.post.agg, tstd <= min.post)
    
    Y.post.agg <- data.frame(ID = unit[no.agg == FALSE],
                             Actual = Y.post.agg$Actual,
                             Synthetic = synth.mat[no.agg == FALSE, 1],
                             Time = c(0:(sum(no.agg==FALSE)-1))) 
    
    Y.pre.agg$Treatment <- 0
    Y.post.agg$Treatment <- 1
    
    Y.pre.agg$ID <- "aggregate"
    Y.post.agg$ID <- "aggregate"
    
    Y.actual <- rbind(Y.pre.agg, Y.post.agg)
    Y.actual$Tdate <- 0
    
    plot.type <- "unit-time"
    I <- 1
    treated.reception <- data.frame(ID="aggregate", Tdate = 1/2)
    toplot <- Y.actual
    toplot$Time <- toplot$Time + 1
    
  } else {
    # Merge synthetic
    names <- strsplit(rownames(synth.mat), "\\.")
    unit <- unlist(lapply(names, "[[", 1))
    period <- unlist(lapply(names, "[[", 2))
    
    synth <- data.frame(ID = unit, Time = period, Synthetic = synth.mat)
    toplot <- merge(Y.actual, synth, by = c("ID", "Time"), all = FALSE) # keep only treated units
  }
 
  return(list(toplot=toplot, treated.reception=treated.reception, plot.type=plot.type)) 
}


