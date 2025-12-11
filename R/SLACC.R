SLACC = function(dat, mod = NULL, L = 5, batch = NULL, include_diag = T, init = NULL,
                  maxIter = 30, eps = 1e-3, U_maxIter = 100, U_eps = 1e-3, 
                  lambda_U = NULL, tau = NULL, lambda_BIC = NULL){
  n = nrow(dat); p = ncol(dat); V = (sqrt(1+8*p)-1)/2;
  
  if (is.null(batch)){ batch = as.factor(rep("group 1", n)) } 
  ni = as.integer(table(batch))
  groups = split(seq_len(n), batch)
  M = length(groups)
  batch_levels = levels(batch) 
  
  if (is.null(mod)){
    if ( M==1 ){ X = model.matrix(~1) } 
    else{ X = cbind(model.matrix(~batch-1)) }
  } else{
    if ( M==1 ){ X = model.matrix(~mod) } 
    else{ X = cbind(model.matrix(~batch-1), model.matrix(~mod-1)) }    
  }
  
  q = ncol(X)
  if (include_diag){ nonzero=1:p }
  else{ nonzero = which(Ltrans(diag(V),d = T)==0) }
  
  p0 = length(nonzero)
  
  if ( is.null(tau) ){ tau = 0.5*sqrt(log(V*L)/n) }
  if ( is.null(lambda_U) ){ lambda_U = log(n*V*L) }
  if ( is.null(lambda_BIC)){ lambda_BIC = log(n)+ 2*gamma*log(p0) }
  
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
  active = 1:L
  la = length(active)
  while (Iter < maxIter){
    Iter=Iter+1
    cat(paste0("Iteration=",Iter), "\n")
    
    #Preliminary
    StS = crossprod(S[nonzero, active ,drop=FALSE])
    Rinv = chol2inv(chol(R[active, active]+1e-8*diag(la)))                    
    Q_list = vector("list", M)
    SigmaAinv_list = vector("list", M)    
    A = matrix(0,n,L)
    for (g in 1:M) {
      idx = groups[[g]]
      d_g = sqrt(pmax(as.numeric(sigma2_g[g,active]), 1e-8)) 
      Dinv_g = diag(1/d_g, la, la)
      SigmaAinv_list[[g]] = matrix(0, L, L)
      SigmaAinv_list[[g]][active,active] = Dinv_g %*% Rinv %*% Dinv_g
      Q_list[[g]] = matrix(0,L,L)
      Q_list[[g]][active,active] = chol2inv(chol(SigmaAinv_list[[g]][active,active] + StS / phi2_g[g]+1e-8*diag(la)))
      A[idx, active] = (X[idx,,drop=FALSE] %*% B[,active,drop=FALSE] %*% SigmaAinv_list[[g]][active,active] + (dat[idx,nonzero,drop=FALSE] %*% S[nonzero,active,drop=FALSE])/phi2_g[g]) %*% Q_list[[g]][active,active]
    }
    
    #M step - update U
    prep = prepare_elements(dat, A=A[,active], U=U[,active], L=la, phi2=phi2_g, tau=tau, ni=ni)
    if (lambda_U>0){
      Q_list_act = lapply(Q_list, function(Qg) Qg[active, active, drop = FALSE])
      Utemp = bilinear_admm(Y = prep$Y, A = prep$X, w = prep$subj_wts, Q = Q_list_act, groups=groups, C = prep$B_wts, U0=U[,active], V0=U[,active], lambda = lambda_U/2, maxit = U_maxIter, tol = U_eps, include_diag=include_diag)$U
      U = matrix(0, V, L)
      U[,active] = Utemp
    } else if (lambda_U==0){
      U = bilinear_als(Y = prep$Y, A = prep$X, w = prep$subj_wts, Q = Q_list, groups = groups, U0 = U, V0 = U, include_diag = include_diag, maxit = U_maxIter, tol = U_eps )$U
    }
    
    active = which(colSums(U^2) != 0)
    la = length(active)
    if (la < L) { A[, -active] = 0 }

    S = foreach(l=1:L,.combine="cbind")%do%{ Ltrans(tcrossprod(U[,l])) }
    
    #M step - update B
    B = matrix(0, nrow = q, ncol = L) 
    if (la > 0) {
      H_act = matrix(0, q * la, q * la)
      b_act = numeric(q * la)
      
      for (g in seq_along(groups)) {
        idx <- groups[[g]]
        Xg = X[idx,,drop = FALSE]          # n_g × q
        Ag = A[idx,,drop = FALSE]          # n_g × L
        
        Ag_act = Ag[, active, drop = FALSE]   # n_g × la
        SigA_g = SigmaAinv_list[[g]][active, active, drop = FALSE]  # la × la
        
        H_act = H_act + kronecker(t(SigA_g), crossprod(Xg))
        b_act = b_act + as.vector(crossprod(Xg, Ag_act %*% SigA_g))
      }
      
      H_act_reg = H_act + 1e-8 * diag(q * la)
      vecB_act = solve(H_act_reg, b_act)
      B[, active] = matrix(vecB_act, nrow = q, ncol = la, byrow = FALSE)
    }
    
    #M step - update sigma2
    R_accum = matrix(0, L, L); 
    for (g in 1:M) {
      idx = groups[[g]]; ng = length(idx)
      RgA = A[idx, , drop=FALSE] - X[idx, , drop=FALSE] %*% B   
      Sigma_hat_g = crossprod(RgA)/ng + Q_list[[g]]             
      sigma2_g[g, active] = diag(Sigma_hat_g)[active]
      sigma2_g[g, -active] = 0
      Dg_act = diag(sqrt(sigma2_g[g, active]), la, la)
      Cg_act = solve(Dg_act, Sigma_hat_g[active,active, drop=FALSE])
      Cg_act = Cg_act %*% solve(Dg_act)
      Cg_act = (Cg_act + t(Cg_act))/2
      Cg = matrix(0, L, L)
      Cg[active, active] = Cg_act
      R_accum = R_accum + ng * Cg
    }
    R_raw = R_accum / n
    R = cov2cor(R_raw+1e-8*diag(L))
    
    if (la < L) { 
      sigma2_g[,-active] = 0
      R[-active, ] = R[,-active] = 0
      diag(R) = 1
      }
    
    #M step - update phi2
    for (g in 1:M) {
      idx = groups[[g]]; ng = length(idx)
      Rg = dat[idx,nonzero,drop=FALSE] - tcrossprod(A[idx,,drop=FALSE], S[nonzero,,drop=FALSE])
      phi2_g[g] = (sum(Rg^2) + ng * sum(diag(S[nonzero,,drop=FALSE] %*% Q_list[[g]] %*% t(S[nonzero,,drop=FALSE])))) / (ng * p0)
    }
    
    order = align_loadings(U=U_prev,U,method = "corr")
    Uhat = order$Uhat_aligned
    
    if ( norm(U_prev-Uhat, type="2")/sqrt(V*L) < eps){ break }
    else{ U_prev=U }
  }
  
  #Post-hoc scaling
  for (l in 1:L){
    sc = sqrt(sum(U[,l]^2))
    if (sc > 0){
      U[,l] = U[,l]/sc
      A[,l] = A[,l]*sc^2
      B[,l] = B[,l]*sc^2
      sigma2_g[,l] = sigma2_g[,l]*sc^4
    }
  }
  S = foreach(l=1:L, .combine="cbind") %do% { Ltrans(tcrossprod(U[,l])) }
  
  resid_mean_list = vector("list", M)
  names(resid_mean_list) = batch_levels
  
  for (g in seq_len(M)) {
    idx = groups[[g]]
    if (length(idx) == 0) next
    
    E_g = dat[idx, nonzero, drop = FALSE] - tcrossprod(A[idx, , drop = FALSE], S[nonzero, , drop = FALSE])
    resid_mean_list[[g]] = colMeans(E_g)
  }
  
  
  ll = logLikSLACC_batch(dat[,nonzero,drop=FALSE], X, B[,active], S[nonzero,active,drop=FALSE], R[active,active], sigma2_by_batch = sigma2_g[,active], phi2_by_batch = phi2_g, batch = batch)
  nparam = sum(U!=0) + sum(B!=0) + sum(sigma2_g!=0) + M + sum(Ltrans(R[active,active], d = F)!=0)

  estimates = list(A = A, S = S, U = U, B = B, R = R, sigma2 = sigma2_g, phi2 = phi2_g, resid_means = resid_mean_list)
  
  input = list(X = X, L = L, batch = batch, lambda_U = lambda_U, lambda_BIC = lambda_BIC, tau = tau, gamma = gamma,
               maxIter = maxIter, eps = eps, U_maxIter = U_maxIter, U_eps = U_eps, include_diag = include_diag)
  
  measure = list(logLik = ll, df = nparam)

  return(list(estimates = estimates, input = input, measure = measure))
}