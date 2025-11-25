harmonize_SLACC_external <- function(dat_new, mod_new, batch_new, fit, eps = 1e-10){
  # ---- pull from fit ----
  est = fit$estimates
  harm = fit$estimates_harmonization
  nonzero = fit$nonzero
  if (is.null(nonzero)) nonzero = seq_len(ncol(dat_new))
  B = est$B                  # q x L
  S = est$S[nonzero, , drop = FALSE]  # p0 x L
  R = est$R
  sigma2_g = est$sigma2          # M x L
  phi2_g = est$phi2            # length M
  batch_levels = fit$batch_levels
  q = nrow(B)
  L = ncol(B)
  M = length(batch_levels)
  
  if (M <= 1 || is.null(harm)) {
    # 배치가 1개면 조화 필요 없음: 원자료 그대로(비활성 엣지는 그대로 유지)
    out <- dat_new
    return(list(dat_harmonized = out, A_new = NULL, A_harmonized = NULL))
  }
  
  # ---- design matrix for new data (훈련 설계와 동일한 열 순서/정의) ----
  batch_new <- factor(batch_new, levels = batch_levels)
  if (any(is.na(batch_new))) stop("새 배치 레벨이 학습時 batch_levels에 없습니다.")
  if (M == 1) {
    X_new = model.matrix(~mod_new)
  } else {
    X_new = cbind(model.matrix(~batch_new - 1),
                   model.matrix(~mod_new - 1))
  }
  if (ncol(X_new) != q) stop("X_new의 열 수가 학습된 B의 행 수(q)와 다릅니다.")
  
  n_new = nrow(dat_new)
  p_new = ncol(dat_new)
  p0 = length(nonzero)
  Y_nn = dat_new[, nonzero, drop = FALSE]
  
  # ---- one-step E-update for A on external data (re-fit 없이) ----
  StS = crossprod(S)                             # L x L
  Rinv = chol2inv(chol(R + eps * diag(L)))
  groups_new = split(seq_len(n_new), batch_new)
  
  A_new = matrix(NA, nrow = n_new, ncol = L)
  for (g in seq_len(M)) {
    idx = groups_new[[g]]
    if (length(idx) == 0) next
    d_g = sqrt(pmax(as.numeric(sigma2_g[g, ]), 1e-12))
    Dinvg = diag(1/d_g, L, L)
    SAinv_g = Dinvg %*% Rinv %*% Dinvg
    Q_g = chol2inv(chol(SAinv_g + StS / phi2_g[g] + eps * diag(L)))
    A_new[idx, ] = ( X_new[idx, , drop = FALSE] %*% B %*% SAinv_g +
                        (Y_nn[idx, , drop = FALSE] %*% S) / phi2_g[g] ) %*% Q_g
  }
  
  # ---- harmonization (훈련과 동일한 규칙, 학습된 타깃 사용) ----
  sigma2_star = harm$sigma2_star   # length L
  phi2_star = harm$phi2_star     # scalar
  center_vec = harm$center_vec    # length L
  
  gamma_hat = B[1:M, , drop = FALSE]                      # M x L
  gamma_cent = sweep(gamma_hat, 2, center_vec, "-")       # M x L
  XB_mean_harm = X_new %*% B - X_new[, 1:M, drop = FALSE] %*% gamma_cent  # n x L
  
  A_harmonized = matrix(NA, n_new, L)
  for (g in seq_len(M)) {
    idx = groups_new[[g]]
    if (length(idx) == 0) next
    A_resid = A_new[idx, , drop = FALSE] - X_new[idx, , drop = FALSE] %*% B
    sf = sqrt(sigma2_star) / sqrt(pmax(as.numeric(sigma2_g[g, ]), 1e-12))
    A_resid_h = sweep(A_resid, 2, sf, "*")
    A_harmonized[idx, ] = XB_mean_harm[idx, ] + A_resid_h
  }
  
  # 관측 잔차 스케일 + 재조합 (nonzero 위치만 채움)
  dat_harm = matrix(0, n_new, p_new)
  for (g in seq_len(M)) {
    idx = groups_new[[g]]
    if (length(idx) == 0) next
    E_g = Y_nn[idx, , drop = FALSE] - tcrossprod(A_new[idx, , drop = FALSE], S)
    
    ##Trial code
    E_g = E_g-rep(1, length(idx))%*%t(colMeans(E_g))
    ####
                  
    
    sphi = sqrt(phi2_star / pmax(as.numeric(phi2_g[g]), 1e-12))
    E_g_h = sphi * E_g
    Signal_h = tcrossprod(A_harmonized[idx, , drop = FALSE], S)
    dat_harm[idx, nonzero] = Signal_h + E_g_h
  }
  
  return(list(
    dat_harmonized = dat_harm,
    A_new = A_new,
    A_harmonized = A_harmonized
  ))
}