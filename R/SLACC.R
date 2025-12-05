SLACC = function(dat, mod = NULL, L = 5, batch = NULL, include_diag = T, init = NULL,
                  maxIter = 30, eps = 1e-3, U_maxIter = 100, U_eps = 1e-3, 
                  lambda_U = NULL, tau = NULL, lambda_BIC = NULL, gamma = 0.5){
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
    if (lambda_U>0){
      U = bilinear_admm(Y = prep$Y, A = prep$X, w = prep$subj_wts, Q = Q_list, groups=groups, C = prep$B_wts, U0=U, V0=U, lambda = lambda_U/2, maxit = U_maxIter, tol = U_eps, include_diag=include_diag)$U
    } else if (lambda_U==0){
      U = bilinear_als(Y = prep$Y, A = prep$X, w = prep$subj_wts, Q = Q_list, groups = groups, U0 = U, V0 = U, include_diag = include_diag, maxit = U_maxIter, tol = U_eps )$U
    }
    
    for (l in 1:L){
      sc = sqrt(sum(U[,l]^2))
      U[,l] = U[,l]/sc; A[,l] = A[,l]*sc^2
    }
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
    
    order = align_loadings(U=U_prev,U,method = "corr")
    Uhat = order$Uhat_aligned
    
    if ( norm(U_prev-Uhat, type="2")/sqrt(V*L) < eps){ break }
    else{ U_prev=U }
  }
  
  ll = logLikSLACC_batch(dat[,nonzero,drop=FALSE], X, B, S[nonzero,,drop=FALSE], R, sigma2_by_batch = sigma2_g, phi2_by_batch = phi2_g, batch = batch)
  nparam = sum(U!=0) + q*L + M*L + M + L*(L-1)/2
  BIC = -2*ll + lambda_BIC*nparam
  
  estimates = list(A = A, S = S, U = U, B = B, R = R, sigma2 = sigma2_g, phi2 = phi2_g)
  
  input = list(X = X, L = L, batch = batch, lambda_U = lambda_U, lambda_BIC = lambda_BIC, tau = tau, gamma = gamma,
               maxIter = maxIter, eps = eps, U_maxIter = U_maxIter, U_eps = U_eps, include_diag = include_diag)
  
  measure = list(logLik = ll, BIC = BIC, df = nparam)
  
  return(list(estimates = estimates, input = input, measure = measure))
}