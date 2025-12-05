SLACC_harmonize = function(dat, batch, fit){
  n = nrow(dat); p = ncol(dat)
  batch = fit$input$batch
  ni = as.integer(table(batch))
  groups = split(seq_len(n), batch)
  M = length(groups)
  A = fit$estimates$A
  S = fit$estimates$S
  B = fit$estimates$B
  sigma2 = fit$estimates$sigma2
  phi2 = fit$estimates$phi2
  X = fit$input$X
  
  #Harmonization step
  estimates_harmonization = NULL
  dat_harmonized = NULL
  dat_harmonized = matrix(0, nrow = n, ncol = p)
  
  wi = numeric(n)
  for (g in seq_len(M)) wi[groups[[g]]] = 1/(2*pmax(as.numeric(phi2[g]), 1e-8))
  wi = wi / sum(wi)
  
  w_g = vapply(seq_len(M), function(g) sum(wi[groups[[g]]]), 0.0)
  w_g = pmax(w_g, 1e-8); w_g = w_g/sum(w_g)
  
  sigma2_star = as.numeric(colSums(sigma2 * w_g))
  phi2_star = sum(w_g * phi2)    
  
  gamma_hat = B[1:M,,drop=FALSE]   
  A_harmonized = matrix(NA, nrow = n, ncol = L)
  
  #Mean harmonization
  center_vec=as.numeric(crossprod(w_g, gamma_hat))  
  gamma_cent = sweep(gamma_hat, 2, center_vec, "-")    
  XB_mean_harmonized = X %*% B - X[,1:M,drop=F] %*% gamma_cent   
  
  #Latent residual harmonization
  for (g in 1:M) {
    idx = groups[[g]]
    A_resid_idx = A[idx,,drop=FALSE] - X[idx,,drop=FALSE] %*% B   
    
    sf = sqrt(sigma2_star)/sqrt(as.numeric(sigma2[g, ])) 
    A_resid_harmonized_idx = sweep(A_resid_idx, 2, sf, FUN = "*")
    
    A_harmonized[idx, ] = XB_mean_harmonized[idx,] + A_resid_harmonized_idx
  }
  
  #Residual harmonization
  # cm_g = matrix(NA, M, p0)
  for (g in 1:M) {
    idx = groups[[g]]
    E_g = dat[idx,nonzero,drop=FALSE] - tcrossprod(A[idx,,drop=FALSE], S[nonzero,,drop=FALSE])
    
    # cm_g[g,] = colMeans(E_g)
    # E_g = E_g-rep(1, length(idx))%*%t(cm_g[g,])
    
    sphi = sqrt(phi2_star/phi2[g])
    E_g_harmonized = sphi * E_g
    dat_harmonized[idx,nonzero] = tcrossprod(A_harmonized[idx,,drop=FALSE], S[nonzero,,drop=FALSE]) + E_g_harmonized
  }
  
  estimates_harmonization=list(A_harmonized = A_harmonized, sigma2_star = sigma2_star, phi2_star = phi2_star, center_vec = center_vec)
  
  return(list(estimates_harmonization = estimates_harmonization, dat_harmonized = dat_harmonized, nonzero=nonzero))
}