SLACC_harmonize = function(dat, fit, mod = NULL, batch = NULL) {
  n = nrow(dat)
  p = ncol(dat)
  
  est = fit$estimates
  A_tr = est$A
  S_hat = est$S
  B_hat = est$B
  U_hat = est$U
  R_hat = est$R
  sigma2_g = est$sigma2
  phi2_g = est$phi2
  L = ncol(B_hat)
  
  inp = fit$input
  X_tr = inp$X
  batch_tr = inp$batch
  include_diag = if (!is.null(inp$include_diag)) inp$include_diag else TRUE
  
  ## --------- nonzero ---------
  if (include_diag) {
    nonzero = seq_len(p)
  } else {
    V = (sqrt(1 + 8*p) - 1) / 2
    nonzero = which(Ltrans(diag(V), d = TRUE) == 0)
  }
  p0 = length(nonzero)
  
  batch_tr = as.factor(batch_tr)
  groups_tr = split(seq_along(batch_tr), batch_tr)
  M = length(groups_tr)
  batch_lvls = levels(batch_tr)
  
  wi_tr = numeric(length(batch_tr))
  for (g in seq_len(M)) {
    wi_tr[groups_tr[[g]]] = 1 / (2 * pmax(as.numeric(phi2_g[g]), 1e-8))
  }
  wi_tr = wi_tr / sum(wi_tr)
  
  ## batch-level weights
  w_g = vapply(seq_len(M), function(g) sum(wi_tr[groups_tr[[g]]]), 0.0)
  w_g = pmax(w_g, 1e-8)
  w_g = w_g / sum(w_g)
  
  ## --------- active set in latent dimension ---------
  active = which(colSums(U_hat^2) != 0)
  la = length(active)
  if (la == 0) {
    stop("All latent factors are zero (no active column in U). Harmonization not well-defined.")
  }
  
  ## --------- construct X_new, batch_new, groups_new ---------
  if (is.null(mod) && is.null(batch)) {
    ## internal harmonization
    if (n != nrow(X_tr)) {
      stop("For internal harmonization, 'dat' must have same n as training data.")
    }
    X_new = X_tr
    batch_new  = batch_tr
  } else if (!is.null(mod) && !is.null(batch)) {
    ## external harmonization
    batch_new = factor(batch, levels = batch_lvls)
    if (any(is.na(batch_new))) {
      stop("External 'batch' has levels not seen in training.")
    }
    if (length(batch_new) != n) {
      stop("Length of 'batch' must match nrow(dat).")
    }
    
    if (M == 1) {
      if (is.null(mod)) {
        X_new = model.matrix(~ 1)
      } else {
        X_new = model.matrix(~ mod)
      }
    } else {
      if (is.null(mod)) {
        X_new = cbind(model.matrix(~ batch_new - 1))
      } else {
        X_new = cbind(model.matrix(~ batch_new - 1),
                      model.matrix(~ mod - 1))
      }
    }
    
    if (ncol(X_new) != nrow(B_hat)) {
      stop("Design matrix for external data has incompatible dimension with training B.")
    }
  } else {
    stop("You must either provide both 'mod' and 'batch' (external mode), or neither (internal mode).")
  }
  
  groups_new = split(seq_len(n), batch_new)
  
  ## --------- E-step: A_new (latent scores for new data) ---------
  A_new = matrix(0, nrow = n, ncol = L)
  
  Rinv = chol2inv(chol(R_hat[active, active, drop = FALSE] + 1e-8 * diag(la)))
  StS  = crossprod(S_hat[nonzero, active, drop = FALSE])
  
  for (g in seq_len(M)) {
    idx = groups_new[[g]]
    if (length(idx) == 0) next
    
    d_g = sqrt(pmax(as.numeric(sigma2_g[g, active]), 1e-8))
    Dinv_g  = diag(1 / d_g, la, la)
    SigmaAinv_act = Dinv_g %*% Rinv %*% Dinv_g               # la x la
    
    Q_g_act = chol2inv(chol(SigmaAinv_act + StS / phi2_g[g] +
                              1e-8 * diag(la)))              # la x la
    
    A_new[idx, active] =
      ( X_new[idx, , drop = FALSE] %*% B_hat[, active, drop = FALSE] %*% SigmaAinv_act +
          (dat[idx, nonzero, drop = FALSE] %*%
             S_hat[nonzero, active, drop = FALSE]) / phi2_g[g]
      ) %*% Q_g_act
  }
  
  ## --------- pooled variance sigma2*, phi2* (active only) ---------
  sigma2_star_act = colSums(sigma2_g[, active, drop = FALSE] * w_g)
  phi2_star = sum(w_g * phi2_g)
  
  ## --------- mean harmonization in latent space ---------
  gamma_hat = B_hat[1:M, active, drop = FALSE]    # M x la
  center_vec = as.numeric(crossprod(w_g, gamma_hat))    # 1 x la
  gamma_cent = sweep(gamma_hat, 2, center_vec, "-")     # M x la
  
  XB_all = X_new %*% B_hat                             # n x L
  XB_mean_harmonized = XB_all
  XB_mean_harmonized[, active] =
    XB_all[, active, drop = FALSE] -
    X_new[, 1:M, drop = FALSE] %*% gamma_cent
  
  ## --------- latent residual harmonization ---------
  A_harmonized = matrix(0, nrow = n, ncol = L)
  
  for (g in seq_len(M)) {
    idx = groups_new[[g]]
    if (length(idx) == 0) next
    
    A_resid_idx_act =A_new[idx,active,drop = FALSE] -
      X_new[idx,,drop = FALSE] %*% B_hat[,active,drop = FALSE]
    
    sf_act = sqrt(sigma2_star_act) / sqrt(pmax(as.numeric(sigma2_g[g, active]), 1e-8))
    
    A_resid_harmonized_act = sweep(A_resid_idx_act, 2, sf_act, "*")
    
    A_harmonized[idx, active] = XB_mean_harmonized[idx, active, drop = FALSE] + A_resid_harmonized_act
  }
  
  ## --------- residual harmonization & reconstruct data ---------
  dat_harmonized = matrix(0, nrow = n, ncol = p)
  
  for (g in seq_len(M)) {
    idx = groups_new[[g]]
    if (length(idx) == 0) next
    
    ## original residuals (latent part fixed to A_new)
    E_g = dat[idx, nonzero, drop = FALSE] -
      tcrossprod(A_new[idx, active, drop = FALSE], S_hat[nonzero, active, drop = FALSE])
    
    sphi = sqrt(phi2_star / phi2_g[g])
    E_g_harmonized = sphi * E_g
    
    Signal_g_harmonized =
      tcrossprod(A_harmonized[idx, active, drop = FALSE],
                 S_hat[nonzero, active, drop = FALSE])
    
    dat_harmonized[idx, nonzero] =
      Signal_g_harmonized + E_g_harmonized
  }
  
  ## --------- output ---------
  estimates_harmonization = list(
    A_harmonized = A_harmonized,
    sigma2_star = sigma2_star_act,
    phi2_star = phi2_star,
    center_vec = center_vec,
    active = active
  )
  
  return(list(estimates_harmonization = estimates_harmonization,
              dat_harmonized = dat_harmonized,
              nonzero = nonzero))
}