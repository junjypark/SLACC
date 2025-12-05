SLACC_harmonize = function(dat, fit, mod = NULL, batch = NULL){
  n = nrow(dat); p = ncol(dat)
  
  est = fit$estimates
  A_tr = est$A 
  S = est$S
  B = est$B
  sigma2 = est$sigma2 
  phi2 = est$phi2 
  L = ncol(B)
  
  inp =fit$input
  X_tr =inp$X 
  batch_tr = inp$batch 
  include_diag = if (!is.null(inp$include_diag)) inp$include_diag else TRUE
  
  if (include_diag){ nonzero = seq_len(p) }
  else {
    V = (sqrt(1 + 8*p) - 1)/2
    nonzero = which(Ltrans(diag(V), d = TRUE) == 0)
  }
  p0 = length(nonzero)
  
  batch_tr = as.factor(batch_tr)
  groups_tr = split(seq_along(batch_tr), batch_tr)
  M = length(groups_tr)
  batch_levels = levels(batch_tr)
  
  if (is.null(mod) && is.null(batch)) {
    ## --------- original data ----------
    if (n != nrow(X_tr)) { stop("For internal harmonization, 'dat' must have same n as training data.") }
    X_new = X_tr; batch_new = batch_tr; A_new = A_tr
    
  } else if (!is.null(mod) && !is.null(batch)) {
    ## --------- external data ----------
    
    batch_new = factor(batch, levels = batch_levels)
    if (any(is.na(batch_new))) { stop("External 'batch' has levels not seen in training.") }
    if (length(batch_new) != n) { stop("Length of 'batch' must match nrow(dat).") }
    
    # X_new
    if (M == 1) {
      if (is.null(mod)) { X_new = model.matrix(~ 1) } 
      else { X_new = model.matrix(~ mod) }
    } else {
      if (is.null(mod)) { X_new = cbind(model.matrix(~ batch_new - 1)) } 
      else { X_new = cbind(model.matrix(~ batch_new - 1), model.matrix(~ mod - 1)) }
    }
    
    if (ncol(X_new) != nrow(B)) {
      stop("Design matrix for external data has incompatible dimension with training B.")
    }
    
    ## ---- E-step ----
    R_mat = est$R
    Rinv = chol2inv(chol(R_mat + 1e-8 * diag(L)))
    StS = crossprod(S[nonzero, , drop = FALSE])
    
    A_new = matrix(NA, nrow = n, ncol = L)
    
    groups_new = split(seq_len(n), batch_new)  # levels(batch_new) == batch_levels
    
    for (g in seq_len(M)) {
      idx = groups_new[[g]]
      if (length(idx) == 0) next
      
      d_g = sqrt(pmax(as.numeric(sigma2[g, ]), 1e-8))
      Dinv_g = diag(1/d_g, L, L)
      SigmaAinv_g = Dinv_g %*% Rinv %*% Dinv_g
      phi2_g = phi2_vec[g]
      
      Q_g = chol2inv(chol(SigmaAinv_g + StS / phi2_g + 1e-8 * diag(L)))
      
      A_new[idx, ] = (X_new[idx,,drop = FALSE] %*% B %*% SigmaAinv_g+
                        (dat[idx,nonzero,drop = FALSE] %*% S[nonzero,,drop = FALSE])/phi2_g )%*%Q_g
    }
  } else { stop("You must either provide both 'mod' and 'batch' (external mode), or neither (internal mode).") }
  
  
  # subject-level weight (training 기준)
  wi_tr = numeric(length(batch_tr))
  for (g in seq_len(M)) {
    wi_tr[groups_tr[[g]]] = 1/(2*pmax(as.numeric(phi2_vec[g]), 1e-8))
  }
  wi_tr = wi_tr / sum(wi_tr)
  
  # batch-level weight
  w_g = vapply(seq_len(M), function(g) sum(wi_tr[groups_tr[[g]]]), 0.0)
  w_g = pmax(w_g, 1e-8); w_g = w_g / sum(w_g)
  
  sigma2_star = as.numeric(colSums(sigma2 * w_g)) 
  phi2_star = sum(w_g * phi2_vec)
  
  # mean harmonization용 gamma_hat (batch effect in latent space)
  gamma_hat = B[1:M, , drop = FALSE]    # M × L (batch dummy coefficient)
  
  ##harmonization: A_new, dat
  if (is.null(mod) && is.null(batch)) {
    batch_new = batch_tr
    X_new = X_tr
    groups_new = split(seq_len(n), batch_new)
  } else {
    groups_new = split(seq_len(n), batch_new)
  }
  
  A_harmonized = matrix(NA, nrow = n, ncol = L)
  
  ## (1) mean harmonization
  center_vec = as.numeric(crossprod(w_g, gamma_hat))       # 1×L
  gamma_cent = sweep(gamma_hat, 2, center_vec, "-")        # M×L
  XB_mean_harmonized = X_new%*%B - X_new[,1:M,drop = FALSE] %*% gamma_cent
  
  ## (2) latent residual harmonization
  for (g in seq_len(M)) {
    idx = groups_new[[g]]
    if (length(idx) == 0) next
    
    A_resid_idx = A_new[idx, , drop = FALSE] - X_new[idx, , drop = FALSE] %*% B
    
    sf = sqrt(sigma2_star) / sqrt(pmax(as.numeric(sigma2[g, ]), 1e-8))
    A_resid_harmonized_idx = sweep(A_resid_idx, 2, sf, "*")
    A_harmonized[idx, ] = XB_mean_harmonized[idx, ] + A_resid_harmonized_idx
  }
  
  ## (3) residual harmonization
  dat_harmonized = matrix(0, nrow = n, ncol = p)
  
  for (g in seq_len(M)) {
    idx = groups_new[[g]]
    if (length(idx) == 0) next
    
    E_g = dat[idx, nonzero, drop = FALSE] - tcrossprod(A_new[idx,,drop = FALSE], S[nonzero,,drop = FALSE])
    sphi = sqrt(phi2_star / phi2_vec[g])
    E_g_harmonized = sphi * E_g
    
    Signal_g_harmonized = tcrossprod(A_harmonized[idx,,drop = FALSE], S[nonzero,,drop = FALSE])
    
    dat_harmonized[idx, nonzero] = Signal_g_harmonized + E_g_harmonized
  }
  
  estimates_harmonization = list(A_harmonized = A_harmonized,sigma2_star = sigma2_star,
                                 phi2_star = phi2_star, center_vec = center_vec
  )
  
  return(list(estimates_harmonization = estimates_harmonization,
              dat_harmonized = dat_harmonized, nonzero = nonzero))
}