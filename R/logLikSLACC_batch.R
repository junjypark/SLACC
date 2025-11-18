logLikSLACC_batch = function(dat, mod, B, S, R, sigma2_by_batch, phi2_by_batch, batch){
  n = nrow(dat); p = ncol(dat); L = ncol(S)
  batch  = as.factor(batch)
  groups = split(seq_len(n), batch)
  M = length(groups)
  
  Rinv = chol2inv(chol(R))
  logdetR = as.numeric(determinant(R, logarithm = TRUE)$modulus)
  
  ll = -(n*p/2) * log(2*pi)
  StS = crossprod(S)
  
  for (g in 1:M) {
    idx = groups[[g]]; ng = length(idx)
    phi2_g = as.numeric(phi2_by_batch[g])
    sig2_g = if (is.matrix(sigma2_by_batch)) sigma2_by_batch[g, ] else sigma2_by_batch[[g]]
    dg_sd = sqrt(pmax(as.numeric(sig2_g), 1e-12))
    Dinv_g = diag(1/dg_sd, L, L)
    SAinv_g = Dinv_g %*% Rinv %*% Dinv_g   
    M_g = SAinv_g + StS/phi2_g
    Rchol = chol(M_g)
    logdet_M_g = 2 * sum(log(diag(Rchol)))
    logdet_SigmaA_g = sum(log(sig2_g)) + logdetR
    logdet_SigY_g = p*log(phi2_g) + logdet_M_g + logdet_SigmaA_g
    
    Rg = dat[idx,,drop=FALSE] - (mod[idx,,drop=FALSE] %*% B) %*% t(S) 
    RtR = crossprod(Rg)                  
    Q_g = chol2inv(chol(M_g))
    quad_g = (1/phi2_g) * sum(Rg^2) - (1/phi2_g^2) * sum( Q_g * ( t(S) %*% RtR %*% S ) )
    
    ll = ll - 0.5 * (ng * logdet_SigY_g + quad_g)
  }
  ll
}