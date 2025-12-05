#' Bilinear ALS update for U
#'
#' @export
bilinear_als = function(Y, A, w,
                          U0 = NULL, V0 = NULL,
                          Q = NULL, groups = NULL,
                          include_diag = TRUE,
                          maxit = 50, tol = 1e-3,
                          gamma = 1e-6) {
  
  n = dim(Y)[1]
  Vdim = dim(Y)[2]
  Ldim = dim(A)[2]
  
  stopifnot(length(dim(Y)) == 3L, dim(Y)[3] == Vdim)
  stopifnot(length(dim(A)) == 3L,
            dim(A)[1] == n,
            dim(A)[2] == Ldim,
            dim(A)[3] == Ldim)
  stopifnot(length(w) == n)
  
  wsum = sum(w)
  
  ## ---- Q -> Qbar ----
  if (is.null(Q)) {
    Qbar = matrix(0, Ldim, Ldim)
  } else if (is.matrix(Q)) {
    stopifnot(all(dim(Q) == c(Ldim, Ldim)))
    Qbar = Q
  } else if (is.list(Q)) {
    stopifnot(!is.null(groups))
    M = length(Q)
    stopifnot(length(groups) == M)
    
    Qbar = matrix(0, Ldim, Ldim)
    for (g in seq_len(M)) {
      idx_g = groups[[g]]
      stopifnot(all(idx_g %in% seq_len(n)))
      wg = sum(w[idx_g])
      Qbar = Qbar + wg * Q[[g]]
    }
    Qbar = Qbar / wsum
  } else {
    stop("Q must be NULL, a matrix, or a list of matrices.")
  }
  
  ## ---- initialize ----
  if (is.null(U0)) {
    Ybar = apply(Y, c(2, 3), mean)
    ev   = eigen((Ybar + t(Ybar)) / 2, symmetric = TRUE)
    U = ev$vectors[, seq_len(Ldim), drop = FALSE]
    V = U
  } else {
    U = U0
    V = if (is.null(V0)) U else V0
  }
  
  ## ---- R array -> C++ cube: (V x V x n), (L x L x n) ----
  #   R: dim(Y) = c(n, V, V)
  #   C++: Y_cube: V x V x n
  Y_cpp = aperm(Y, c(2, 3, 1))
  A_cpp = aperm(A, c(2, 3, 1))
  
  res = bilinear_als3_cpp(
    Y_cpp, A_cpp, w,
    U, V,
    Qbar,
    include_diag,
    maxit,
    tol,
    gamma
  )
  
  list(U = res$U, V = res$V)
}