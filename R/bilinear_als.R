bilinear_als3 = function(Y, A, w,
                         U0 = NULL, V0 = NULL,
                         Q = NULL, groups = NULL,
                         include_diag = TRUE,
                         maxit = 50, tol = 1e-3,
                         gamma = 1e-6) {
  # Y: n × V × V  (array)
  # A: n × L × L  (array of A_ij)
  # w: length n   (subject weights, e.g. 1/(2 phi_i^2))
  
  n    = dim(Y)[1]
  Vdim = dim(Y)[2]
  Ldim = dim(A)[2]
  
  stopifnot(length(dim(Y)) == 3, dim(Y)[3] == Vdim)
  stopifnot(length(dim(A)) == 3, dim(A)[1] == n, dim(A)[2] == Ldim, dim(A)[3] == Ldim)
  stopifnot(length(w) == n)
  
  wsum = sum(w)
  
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
  
  ## ---- 초기값 ----
  if (is.null(U0)) {
    Ybar = apply(Y, c(2, 3), mean)
    ev = eigen((Ybar + t(Ybar)) / 2, symmetric = TRUE)
    U = ev$vectors[, seq_len(Ldim), drop = FALSE]
    V = U
  } else {
    U = U0
    V = if (is.null(V0)) U else V0
  }
  
  fro_sq = function(M) sum(M * M)
  
  for (it in seq_len(maxit)) {
    
    ## ---------- U-step ----------
    VtV = crossprod(V)
    GV  = matrix(0, Ldim, Ldim)
    HV  = matrix(0, Vdim, Ldim)
    
    for (i in seq_len(n)) {
      Ai = matrix(A[i, , ], Ldim, Ldim)
      Yi = matrix(Y[i, , ], Vdim, Vdim)
      
      if (!include_diag) {
        a = diag(Ai)           # a_ij (L×1)
        diag_fit = (U * V) %*% a      # diag(U diag(a) V^T) = (U*V) %*% a
        diag(Yi) = as.vector(diag_fit)
      }
      
      GV = GV + w[i] * (Ai %*% VtV %*% t(Ai))
      HV = HV + w[i] * (Yi %*% V %*% t(Ai))
    }
    
    K = GV + wsum * (VtV * Qbar) + gamma * diag(Ldim)
    K = (K + t(K)) / 2
    cholK = chol(K)
    
    U_new = t(backsolve(cholK, forwardsolve(t(cholK), t(HV))))
    
    ## ---------- V-step ----------
    UtU = crossprod(U_new)
    GU  = matrix(0, Ldim, Ldim)
    HU  = matrix(0, Vdim, Ldim)
    
    for (i in seq_len(n)) {
      Ai = matrix(A[i, , ], Ldim, Ldim)
      Yi = matrix(Y[i, , ], Vdim, Vdim)
      
      if (!include_diag) {
        a        = diag(Ai)
        diag_fit = (U_new * V) %*% a   # 여기서도 현재 U_new, V 반영
        diag(Yi) = as.vector(diag_fit)
      }
      
      GU = GU + w[i] * (Ai %*% UtU %*% t(Ai))
      HU = HU + w[i] * (t(Yi) %*% U_new %*% t(Ai))
    }
    
    K2 = GU + wsum * (UtU * Qbar) + gamma * diag(Ldim)
    K2 = (K2 + t(K2)) / 2
    cholK2 = chol(K2)
    
    V_new = t(backsolve(cholK2, forwardsolve(t(cholK2), t(HU))))
    
    r_norm = sqrt(fro_sq(U_new - U) + fro_sq(V_new - V))
    U = U_new
    V = V_new
    
    if (r_norm < tol) break
  }
  
  U_final = (U + V) / 2
  list(U = U_final, V = U_final)
}