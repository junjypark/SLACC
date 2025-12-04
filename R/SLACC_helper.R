Ltrans = function(X,d = TRUE){ X[upper.tri(X,d)]  }

Ltrinv = function(x,V,d = TRUE){ Y = matrix(0,ncol = V,nrow = V);
Y[upper.tri(Y,d)]=x;return(Y + t(Y) - d*diag(diag(Y)))  }

estim_phi2 = function(input, batch, nonzero){
  phi2 = vector()
  ni = table(batch)
  names = unique(batch)
  n.batch = length(ni)
  for (i in 1:n.batch){
    index = which(batch==names[i])
    phi2[i] = var(c(input[index,nonzero,drop=FALSE])) #estim_sigma(input[index,nonzero,drop=FALSE],method="MAD")^2
  }
  return(phi2)
}

HOSVD_initial = function(Y, L, X, batch, nonzero){
  n = nrow(Y)
  p = ncol(Y)
  V = (sqrt(1 + 8*p) - 1) / 2
  
  groups = split(seq_len(n), batch)
  M = length(groups)
  
  ## Y를 (n × V × V) 텐서로 변환
  Yarray = array(NA, dim = c(n, V, V))
  for (i in 1:n){
    Yarray[i,,] = Ltrinv(Y[i,], V)
  }
  Ytensor = rTensor::as.tensor(Yarray)
  
  ## HOSVD는 rank L0 = min(L, V-1)까지만
  L0 = min(L, V-1)
  fit_hosvd = rTensor::hosvd(Ytensor, ranks = c(n, L0, L0))
  
  ## 기본 U, S (L0까지)
  U_ini0 = fit_hosvd$U[[2]][, 1:L0, drop = FALSE]   # V × L0
  U_ini0 = U_ini0 / norm(U_ini0, type = "2")
  
  S_ini0 = foreach(l = 1:L0, .combine = "cbind") %do% {
    Ltrans(tcrossprod(U_ini0[, l]))
  }                                                # p × L0
  
  ## ---- L <= V 인 경우: 여기서 바로 ridge LS 후 리턴 ----
  if (L <= V) {
    S_ini = S_ini0
    U_ini = U_ini0
  } else {
    ## ---- L > V 인 경우: 여분 축(L - L0)을 추가 ----
    L_extra = L - L0
    
    ## (1) U: 기존 U_ini0 + 작은 랜덤 extra 축
    U_extra = matrix(rnorm(V * L_extra, sd = 0.01), nrow = V, ncol = L_extra)
    for (l in 1:L_extra) {
      sc = sqrt(sum(U_extra[, l]^2))
      if (sc > 0) U_extra[, l] = U_extra[, l] / sc
    }
    U_ini = cbind(U_ini0, U_extra)           # V × L
    
    ## (2) S: 기존 S_ini0 + U_extra로부터 만든 extra S
    S_extra = foreach(l = 1:L_extra, .combine = "cbind") %do% {
      Ltrans(tcrossprod(U_extra[, l]))
    }                                        # p × L_extra
    S_ini = cbind(S_ini0, S_extra)           # p × L
  }
  
  ## ---- 공통: 확장된 S_ini(U_ini)에 대해 한 번 더 ridge LS로 A, B, sigma2, R, phi2 추정 ----
  Xmat  <- S_ini[nonzero, , drop = FALSE]         # p0 × L
  XtX   <- crossprod(Xmat)                        # L × L
  eps   <- 1e-6 * mean(diag(XtX))
  XtXr  <- XtX + eps * diag(ncol(XtX))
  
  A_ini <- Y[, nonzero, drop = FALSE] %*% Xmat %*% solve(XtXr)   # n × L
  
  fit_lm = lm(A_ini ~ X - 1)
  B      = coef(fit_lm)                     # q × L
  resid  = residuals(fit_lm)                # n × L
  
  sigma2_ini = foreach(g = 1:M, .combine = "rbind") %do% {
    apply(resid[groups[[g]], , drop = FALSE], 2, var)
  }                                         # M × L
  
  R_ini = cor(resid)                        # L × L
  
  E = Y - A_ini %*% t(S_ini)                # n × p
  phi2_ini = estim_phi2(E, batch, nonzero)  # 길이 M 또는 스칼라
  
  return(list(
    U      = U_ini,
    S      = S_ini,
    A      = A_ini,
    phi2   = phi2_ini,
    B      = B,
    sigma2 = sigma2_ini,
    R      = R_ini
  ))
}