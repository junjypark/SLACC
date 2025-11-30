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
    phi2[i] = estim_sigma(input[index,nonzero,drop=FALSE],method="MAD")^2
  }
  return(phi2)
}

HOSVD_initial=function(Y, L, X, batch, nonzero){
  n = nrow(Y)
  p = ncol(Y)
  V = (sqrt(1+8*p)-1)/2
  
  groups = split(seq_len(n), batch)
  
  Yarray = array(NA, dim = c(n,V,V))
  for (i in 1:n){ Yarray[i,,] = Ltrinv(Y[i,], V) }
  
  Yarray = as.tensor(Yarray)
  
  tmp = tempfile()
  fit_hosvd = rTensor::hosvd(Yarray, ranks = c(n,L,L)) 
  
  U_ini = fit_hosvd$U[[2]]
  U_ini = apply(U_ini, 2, scale)
  S_ini = foreach(l=1:L, .combine="cbind")%do%{ Ltrans(tcrossprod(U_ini[,l])) }
  A_ini = Y %*% S_ini %*% solve(crossprod(S_ini))
  B = lm(A_ini ~ X-1)$coef
  resid=lm(A_ini ~ X-1)$resid
  sigma2_ini=foreach(g=1:length(unique(batch)),.combine="rbind")%do%{
    apply(resid[groups[[g]],], 2,var)
  }
  # R_ini = cor(resid)
  phi2_ini = estim_phi2(Y - A_ini %*% t(S_ini), batch, nonzero)
  # return(list(U = U_ini, S = S_ini, A = A_ini, phi2 = phi2_ini, B = B, sigma2=sigma2_ini, R=R_ini))
  return(list(U = U_ini, S = S_ini, A = A_ini, phi2 = phi2_ini, B = B, sigma2=sigma2_ini))
}

