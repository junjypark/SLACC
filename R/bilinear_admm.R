bilinear_admm = function(Y, A, w, C, lambda,
                               U0 = NULL, V0 = NULL, Q = NULL,
                               rho = 1.0, eta = 1.0, gamma = 1e-6,
                               maxit = 100, tol = 1e-3,
                               snap_to_Z_every = 1) {
  n = dim(Y)[1]; Vdim = dim(Y)[2]; Ldim = dim(A)[2]
  stopifnot(length(dim(Y))==3, dim(Y)[3]==Vdim)
  stopifnot(length(dim(A))==3, dim(A)[1]==n, dim(A)[2]==Ldim, dim(A)[3]==Ldim)
  stopifnot(length(w)==n, all(dim(C)==c(Vdim, Ldim)))
  if (is.null(Q)) Q = matrix(0, Ldim, Ldim)

  # initalize
  if (is.null(U0)) {
    Ybar = apply(Y, c(2,3), mean)
    ev = eigen((Ybar + t(Ybar))/2, symmetric = TRUE)
    U = ev$vectors[, seq_len(Ldim), drop = FALSE]
  } else {
    U = U0
    V = if (is.null(V0)) U else V0
  }

  ZU = U; ZV = V
  WU = matrix(0, Vdim, Ldim); WV = matrix(0, Vdim, Ldim)
  Yeq = matrix(0, Vdim, Ldim)
  wsum = sum(w)
  fro_sq = function(M) sum(M*M)

  for (it in 1:maxit) {
    ## ----- U-step (closed form)
    VtV = crossprod(V)
    GV = matrix(0, Ldim, Ldim); HV = matrix(0, Vdim, Ldim)
    for (i in 1:n) {
      Ai = matrix(A[i,,], Ldim, Ldim)
      Yi = matrix(Y[i,,], Vdim, Vdim)
      GV = GV + w[i] * (Ai %*% VtV %*% t(Ai))
      HV = HV + w[i] * (Yi %*% V %*% t(Ai))
    }
    # K = GV + (rho + eta + gamma) * diag(Ldim) + (wsum/2) * Q
    K = GV + (rho + eta + gamma) * diag(Ldim) + wsum * (VtV*Q)
    K = (K + t(K))/2
    cholK = chol(K)
    T_U = rho*(ZU - WU) + eta*(V - Yeq)
    RHS = HV + T_U
    U = t(backsolve(cholK, forwardsolve(t(cholK), t(RHS))))

    # ZU soft-threshold + dual
    ZU_old = ZU
    ZU = sign(U + WU) * pmax(abs(U + WU) - (lambda/rho) * C, 0)
    WU = WU + U - ZU

    ## ----- V-step (closed form)
    UtU = crossprod(U)
    GU = matrix(0, Ldim, Ldim); HU = matrix(0, Vdim, Ldim)
    for (i in 1:n) {
      Ai = matrix(A[i,,], Ldim, Ldim)
      Yi = matrix(Y[i,,], Vdim, Vdim)
      GU = GU + w[i] * (Ai %*% UtU %*% t(Ai))
      HU = HU + w[i] * (t(Yi) %*% U %*% t(Ai))
    }
    # K = GU + (rho + eta + gamma) * diag(Ldim) + (wsum/2) * Q
    K = GU + (rho + eta + gamma) * diag(Ldim) + wsum * (UtU*Q)
    K = (K + t(K))/2
    cholK = chol(K)
    T_V = rho*(ZV - WV) + eta*(U + Yeq)
    RHS = HU + T_V
    V = t(backsolve(cholK, forwardsolve(t(cholK), t(RHS))))

    # ZV soft-threshold + dual
    ZV_old = ZV
    ZV = sign(V + WV) * pmax(abs(V + WV) - (lambda/rho) * C, 0)
    WV = WV + V - ZV

    # equality dual (U=V)
    Yeq = Yeq + (U - V)

    # snap
    if (snap_to_Z_every > 0 && it %% snap_to_Z_every == 0) { U = ZU; V = ZV }

    # stop
    r_norm = sqrt(fro_sq(U - ZU) + fro_sq(V - ZV) + fro_sq(U - V))
    s_norm = rho * sqrt(fro_sq(ZU - ZU_old) + fro_sq(ZV - ZV_old))
    if (max(r_norm, s_norm) < tol) break
  }
  list(U = ZU, V = ZV, rho = rho)
}