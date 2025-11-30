SLACC2 = function(dat, mod = NULL, L = 5, batch = NULL, maxIter = 20, eps = 1e-3, ADMM_maxIter = 100, ADMM_eps = 1e-3, init = NULL,
                 lambda = NULL, tau = NULL, harmonize = TRUE){
  n = nrow(dat); p = ncol(dat); V = (sqrt(1+8*p)-1)/2;
  
  if (is.null(batch)){ batch = as.factor(rep("group 1", n)) } 
  ni = as.integer(table(batch))
  groups = split(seq_len(n), batch)
  M = length(groups)
  
  if (is.null(mod)){
    if ( M==1 ){ X = model.matrix(~1) } 
    else{ X = cbind(model.matrix(~batch-1)) }
  } else{
    if ( M==1 ){ X = model.matrix(~mod) } 
    else{ X = cbind(model.matrix(~batch-1), model.matrix(~mod-1)) }    
  }

  q = ncol(X)
  nonzero = which(apply(dat,2,sum)!=0)
  p0 = length(nonzero)
  
  if ( is.null(tau) ){ tau = 0.3*sqrt(log(V*L)/n) }
  if ( is.null(lambda) ){ lambda = log(n*p0) }

  #Initialize
  if (is.null(init)){ init = HOSVD_initial(dat, L, X, batch, nonzero) }
  A = init$A; B = init$B; U = U_prev = init$U; S = init$S;
  R= init$R
  phi2_g = init$phi2; sigma2_g = init$sigma2

  
  phi2_g = if (length(phi2_g) == 1) rep(as.numeric(phi2_g), M) else as.numeric(phi2_g)
  if (!is.matrix(sigma2_g) || nrow(sigma2_g) != M || ncol(sigma2_g) != L) {
    sigma2_g = matrix(rep(as.numeric(sigma2_g), length.out = M*L),nrow = M, ncol = L, byrow = TRUE)
  }

  Iter = 0
  while (Iter < maxIter){
    Iter=Iter+1
    cat(paste0("Iteration=",Iter), "\n")
    
    #Preliminary
    StS = crossprod(S[nonzero,,drop=FALSE])
    Rinv = chol2inv(chol(R+1e-8*diag(L)))                    
    Q_list = vector("list", M)
    SigmaAinv_list = vector("list", M)          
    for (g in 1:M) {
      idx = groups[[g]]
      d_g = sqrt(pmax(as.numeric(sigma2_g[g, ]), 1e-8)) 
      Dinv_g = diag(1/d_g, L, L)
      SigmaAinv_list[[g]] = Dinv_g %*% Rinv %*% Dinv_g
      Q_list[[g]] = chol2inv(chol(SigmaAinv_list[[g]] + StS / phi2_g[g]+1e-8*diag(L)))
      A[idx, ] = (X[idx,,drop=FALSE] %*% B %*% SigmaAinv_list[[g]] + (dat[idx,nonzero,drop=FALSE] %*% S[nonzero,,drop=FALSE])/phi2_g[g]) %*% Q_list[[g]]
    }

    #M step - update U
    prep = prepare_elements(dat, A=A, U=U, L=L, phi2=phi2_g, tau=tau, ni=ni)
    w = prep$subj_wts
    w = w/sum(w) #new
    prep$subj_wts = w
    prep$B_wts = (abs(U) <= tau) / tau
    w_g = vapply(groups, function(idx) sum(w[idx]), 0.0)
    Qbar = Reduce("+", Map(function(Qg, wg) wg * Qg, Q_list, as.list(w_g)))
    # Qbar = Qbar / sum(w_g) #new
    U = bilinear_admm(Y = prep$Y, A = prep$X, w = prep$subj_wts, Q = Qbar, C = prep$B_wts, U0=U, lambda = lambda/2, maxit = ADMM_maxIter, tol = ADMM_eps)$U
    S = foreach(l=1:L,.combine="cbind")%do%{ Ltrans(tcrossprod(U[,l])) }

    #M step - update B
    H = matrix(0, q*L, q*L)
    b = numeric(q*L)
    for (g in seq_along(groups)) {
      idx = groups[[g]]
      Xg = X[idx,,drop=FALSE]         # n_g × q
      Ag = A[idx,,drop=FALSE]          # n_g × L
      H = H + kronecker(t(SigmaAinv_list[[g]]), crossprod(Xg))      
      b = b + as.vector(crossprod(Xg, Ag %*% SigmaAinv_list[[g]]))
    }
    vecB = solve(H+1e-8*diag(q*L), b)                      
    B = matrix(vecB, nrow = q, ncol = L, byrow = FALSE)

    #M step - update sigma2
    R_accum = matrix(0, L, L); 
    for (g in 1:M) {
      idx = groups[[g]]; ng = length(idx)
      RgA = A[idx, , drop=FALSE] - X[idx, , drop=FALSE] %*% B   
      Sigma_hat_g = crossprod(RgA)/ng + Q_list[[g]]             
      sigma2_g[g, ] = diag(Sigma_hat_g) 
      Dg = diag(sqrt(sigma2_g[g, ]), L, L)
      Cg = solve(Dg, Sigma_hat_g)
      Cg = Cg %*% solve(Dg)
      Cg = (Cg + t(Cg))/2
      R_accum = R_accum + ng * Cg
    }
    R_raw = R_accum / n
    R = cov2cor(R_raw+1e-8*diag(L))

    #M step - update phi2
    for (g in 1:M) {
      idx = groups[[g]]; ng = length(idx)
      Rg = dat[idx,nonzero,drop=FALSE] - tcrossprod(A[idx,,drop=FALSE], S[nonzero,,drop=FALSE])
      phi2_g[g] = (sum(Rg^2) + ng * sum(diag(S[nonzero,,drop=FALSE] %*% Q_list[[g]] %*% t(S[nonzero,,drop=FALSE])))) / (ng * p0)
    }

    if ( norm(U_prev-U, type="2") < eps){ break }
    else{ U_prev=U }
  }

  ll = logLikSLACC_batch(dat[,nonzero,drop=FALSE], X, B, S[nonzero,,drop=FALSE], R, sigma2_by_batch = sigma2_g, phi2_by_batch = phi2_g, batch = batch)
  nparam = sum(U!=0) + q*L + M*L + M + L*(L-1)/2
  BIC = -2*ll + lambda*nparam

  estimates = list(A = A, S = S, U = U, B = B, R = R, sigma2 = sigma2_g, phi2 = phi2_g)
  
  #Harmonization step
  estimates_harmonization = NULL
  dat_harmonized = NULL
  if (M>1 & harmonize){
    dat_harmonized = matrix(0, nrow = n, ncol = p)
    
    wi = numeric(n)
    for (g in seq_len(M)) wi[groups[[g]]] = 1/(2*pmax(as.numeric(phi2_g[g]), 1e-8))
    wi = wi / sum(wi)
    
    w_g = vapply(seq_len(M), function(g) sum(wi[groups[[g]]]), 0.0)
    w_g = pmax(w_g, 1e-8); w_g = w_g/sum(w_g)
    
    sigma2_star = as.numeric(colSums(sigma2_g * w_g))
    phi2_star = sum(w_g * phi2_g)    
    
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

      sf = sqrt(sigma2_star)/sqrt(as.numeric(sigma2_g[g, ])) 
      A_resid_harmonized_idx = sweep(A_resid_idx, 2, sf, FUN = "*")
      
      A_harmonized[idx, ] = XB_mean_harmonized[idx,] + A_resid_harmonized_idx
    }
    
    #Residual harmonization
    cm_g = matrix(NA, M, p0)
    for (g in 1:M) {
      idx = groups[[g]]
      E_g = dat[idx,nonzero,drop=FALSE] - tcrossprod(A[idx,,drop=FALSE], S[nonzero,,drop=FALSE])
      
      cm_g[g,] = colMeans(E_g)
      E_g = E_g-rep(1, length(idx))%*%t(cm_g[g,])

      sphi = sqrt(phi2_star/phi2_g[g])
      E_g_harmonized = sphi * E_g
      dat_harmonized[idx,nonzero] = tcrossprod(A_harmonized[idx,,drop=FALSE], S[nonzero,,drop=FALSE]) + E_g_harmonized
    }
    
    estimates_harmonization=list(A_harmonized = A_harmonized, cm_g = cm_g, sigma2_star = sigma2_star, phi2_star = phi2_star, center_vec = center_vec)
  } 
  
  return(list(estimates = estimates, estimates_harmonization = estimates_harmonization, logLik = ll, BIC = BIC, dat_harmonized = dat_harmonized, X=X, batch_levels=levels(batch), nonzero=nonzero))
}