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
  
  ## 1) Y를 (n × V × V) 텐서로 복원
  Yarray = array(NA, dim = c(n, V, V))
  for (i in 1:n){
    Yarray[i,,] = Ltrinv(Y[i,], V)
  }
  Ytensor = as.tensor(Yarray)
  
  ## 2) HOSVD로 초기 U 얻기
  fit_hosvd = rTensor::hosvd(Ytensor, ranks = c(n, L, L))
  U_ini = fit_hosvd$U[[2]]
  U_ini = U_ini / norm(U_ini, type = "2")   # 전체 스케일 정규화 (원래 코드 유지)
  
  ## 3) U_ini로부터 S_ini 구성
  S_ini = foreach(l = 1:L, .combine = "cbind") %do% {
    Ltrans(tcrossprod(U_ini[, l]))
  }
  
  ## 4) ridge 붙인 LS로 A_ini 초기값 계산
  Xmat  = S_ini[nonzero, , drop = FALSE]         # p0 × L
  XtX   = crossprod(Xmat)                        # L × L
  eps   = 1e-6 * mean(diag(XtX))                 # 스케일 맞춘 작은 ridge
  XtXr  = XtX + eps * diag(ncol(XtX))            # ridge 추가된 Gram matrix
  
  A_ini = Y[, nonzero, drop = FALSE] %*% Xmat %*% solve(XtXr)  # n × L
  
  ## 5) B, resid, sigma2, R, phi2 초기값
  fit_lm = lm(A_ini ~ X - 1)
  B      = coef(fit_lm)
  resid  = residuals(fit_lm)
  
  sigma2_ini = foreach(g = seq_along(groups), .combine = "rbind") %do% {
    apply(resid[groups[[g]], , drop = FALSE], 2, var)
  }
  
  R_ini = cor(resid)
  
  E = Y - A_ini %*% t(S_ini)
  phi2_ini = estim_phi2(E, batch, nonzero)
  
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