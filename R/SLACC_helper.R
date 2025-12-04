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

# HOSVD_initial=function(Y, L, X, batch, nonzero){
#   n = nrow(Y)
#   p = ncol(Y)
#   V = (sqrt(1+8*p)-1)/2
#   
#   groups = split(seq_len(n), batch)
#   
#   Yarray = array(NA, dim = c(n,V,V))
#   for (i in 1:n){ Yarray[i,,] = Ltrinv(Y[i,], V) }
#   
#   Yarray = as.tensor(Yarray)
#   
#   tmp = tempfile()
#   fit_hosvd = rTensor::hosvd(Yarray, ranks = c(n,L,L)) 
#   
#   U_ini = fit_hosvd$U[[2]]
#   U_ini = U_ini/norm(U_ini, type="2")
#   S_ini = foreach(l=1:L, .combine="cbind")%do%{ Ltrans(tcrossprod(U_ini[,l])) }
#   A_ini = Y[,nonzero] %*% S_ini[nonzero,] %*% solve(crossprod(S_ini[nonzero,]))
#   B = lm(A_ini ~ X-1)$coef
#   resid=lm(A_ini ~ X-1)$resid
#   sigma2_ini=foreach(g=1:length(unique(batch)),.combine="rbind")%do%{
#     apply(resid[groups[[g]],], 2,var)
#   }
#   R_ini = cor(resid)
#   E = Y - A_ini %*% t(S_ini)
#   phi2_ini = estim_phi2(E, batch, nonzero)
#   return(list(U = U_ini, S = S_ini, A = A_ini, phi2 = phi2_ini, B = B, sigma2=sigma2_ini, R=R_ini))
# }

HOSVD_initial = function(Y, L, X, batch, nonzero){
  n = nrow(Y)
  p = ncol(Y)
  V = (sqrt(1 + 8*p) - 1) / 2
  
  groups = split(seq_len(n), batch)
  M = length(groups)
  
  Yarray = array(NA, dim = c(n, V, V))
  for (i in 1:n){
    Yarray[i,,] = Ltrinv(Y[i,], V)
  }
  Yarray = rTensor::as.tensor(Yarray)
  
  L0 = min(L, V)
  fit_hosvd = rTensor::hosvd(Yarray, ranks = c(n, L0, L0))
  
  ## U, S, A, B, sigma2, R, phi2를 먼저 L0에 대해 계산
  U_ini0 = fit_hosvd$U[[2]][, 1:L0, drop = FALSE]   # V × L0
  U_ini0 = U_ini0 / norm(U_ini0, type = "2")
  
  S_ini0 = foreach(l = 1:L0, .combine = "cbind") %do% {
    Ltrans(tcrossprod(U_ini0[, l]))
  }                                                # p × L0
  
  A_ini0 = Y[, nonzero, drop = FALSE] %*% 
    S_ini0[nonzero, , drop = FALSE] %*% 
    solve(crossprod(S_ini0[nonzero, , drop = FALSE]))      # n × L0
  
  fit_lm0 = lm(A_ini0 ~ X - 1)
  B0    = coef(fit_lm0)                     # q × L0
  resid0 = residuals(fit_lm0)               # n × L0
  
  sigma2_ini0 = foreach(g = 1:M, .combine = "rbind") %do% {
    apply(resid0[groups[[g]], , drop = FALSE], 2, var)
  }                                         # M × L0
  
  R_ini0 = cor(resid0)                      # L0 × L0
  
  E0 = Y - A_ini0 %*% t(S_ini0)             # n × p
  phi2_ini = estim_phi2(E0, batch, nonzero) # 길이 M (또는 스칼라)
  
  ## ---- L ≤ V 인 경우: 기존과 동일하게 리턴 ----
  if (L <= V) {
    return(list(
      U      = U_ini0,
      S      = S_ini0,
      A      = A_ini0,
      phi2   = phi2_ini,
      B      = B0,
      sigma2 = sigma2_ini0,
      R      = R_ini0
    ))
  }
  
  L_extra = L - L0
  
  U_extra = matrix(rnorm(V * L_extra, sd = 0.01), nrow = V, ncol = L_extra)
  for (l in 1:L_extra) {
    sc = sqrt(sum(U_extra[, l]^2))
    if (sc > 0) U_extra[, l] = U_extra[, l] / sc
  }
  U_ini = cbind(U_ini0, U_extra)           # V × L
  
  S_extra = foreach(l = 1:L_extra, .combine = "cbind") %do% {
    Ltrans(tcrossprod(U_extra[, l]))
  }                                        # p × L_extra
  S_ini = cbind(S_ini0, S_extra)           # p × L
  
  A_extra = matrix(0, nrow = n, ncol = L_extra)
  A_ini   = cbind(A_ini0, A_extra)         # n × L
  
  B_extra = matrix(0, nrow = ncol(X), ncol = L_extra)
  B       = cbind(B0, B_extra)             # q × L
  
  sigma2_mean_batch = rowMeans(sigma2_ini0)        # 길이 M
  sigma2_extra = matrix(sigma2_mean_batch, nrow = M, ncol = L_extra)
  sigma2_ini   = cbind(sigma2_ini0, sigma2_extra)  # M × L
  
  R_ini = diag(L)
  R_ini[1:L0, 1:L0] = R_ini0
  
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