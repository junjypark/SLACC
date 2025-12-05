# Q_list -> Qbar
compute_Qbar = function(Q, w, groups, L) {
  wsum = sum(w)
  if (is.null(Q)) {
    Qbar = matrix(0, L, L)
  } else if (is.matrix(Q)) {
    stopifnot(all(dim(Q) == c(L, L)))
    Qbar = Q
  } else if (is.list(Q)) {
    stopifnot(!is.null(groups))
    M = length(Q); stopifnot(length(groups) == M)
    Qbar = matrix(0, L, L)
    for (g in seq_len(M)) {
      idx_g = groups[[g]]
      wg    = sum(w[idx_g])
      Qbar  = Qbar + wg * Q[[g]]
    }
    Qbar = Qbar / wsum
  } else {
    stop("Q must be NULL, a matrix, or a list.")
  }
  Qbar
}

bilinear_admm = function(Y, A, w, C, lambda,
                          U0 = NULL, V0 = NULL, Q = NULL,
                          groups = NULL,
                          include_diag = TRUE,
                          rho = 1.0, eta = 1.0, gamma = 1e-6,
                          maxit = 100, tol = 1e-3,
                          snap_to_Z_every = 1) {
  # Y: n x V x V 
  # A: n x L x L
  dimY = dim(Y); n = dimY[1]; V = dimY[2]; stopifnot(dimY[3] == V)
  dimA = dim(A); stopifnot(all(dimA == c(n, dimA[2], dimA[2])))
  L = dimA[2]
  
  stopifnot(length(w) == n, all(dim(C) == c(V, L)))
  
  if (is.null(U0)) {
    Ybar = apply(Y, c(2,3), mean)
    ev   = eigen((Ybar + t(Ybar))/2, symmetric = TRUE)
    U0   = ev$vectors[, seq_len(L), drop = FALSE]
    V0   = U0
  } else {
    if (is.null(V0)) V0 = U0
  }
  
  # Qbar 
  Qbar = compute_Qbar(Q, w, groups, L)
  
  # R array(n,V,V) -> cube(V,V,n)
  Ycube = aperm(Y, c(2,3,1))
  Acube = aperm(A, c(2,3,1))
  
  out = bilinear_admm3_cpp(
    Ycube, Acube,
    w = w,
    C = C,
    lambda = lambda,
    U = U0,
    V = V0,
    Qbar = Qbar,
    include_diag = include_diag,
    rho = rho, eta = eta, gamma_par = gamma,
    maxit = maxit, tol = tol,
    snap_to_Z_every = snap_to_Z_every
  )
  
  out
}