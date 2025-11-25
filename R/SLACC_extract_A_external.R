# dat_new : n_new x p  (외부 데이터, 이미 harmonized 라고 가정)
# X_new   : n_new x q_mod (훈련 때 mod-1 과 동일한 컬럼 구성; 배치 더미 제외)
# fit     : SLACC 학습 결과 객체 (fit$estimates_harmonization 필수)

SLACC_extract_A_external <- function(dat_new, X_new, fit, eps = 1e-10) {
  if (is.null(fit$estimates_harmonization))
    stop("fit$estimates_harmonization 이 필요합니다 (sigma2_star, phi2_star, center_vec).")
  
  # --- pull trained params ---
  est   <- fit$estimates
  S_all <- est$S
  R     <- est$R
  
  hz    <- fit$estimates_harmonization
  sigma2_star <- as.numeric(hz$sigma2_star)
  phi2_star   <- as.numeric(hz$phi2_star)
  
  # 활성 엣지 인덱스
  nonzero <- fit$nonzero
  if (is.null(nonzero)) stop("fit$nonzero 가 없습니다.")
  
  # S (활성 엣지 부분만)
  S <- S_all[nonzero, , drop = FALSE]
  L <- ncol(S)
  
  # B에서 '배치 더미'를 제외한 '생물학적 공변량' 부분만 사용
  B <- est$B
  M <- if (!is.null(fit$batch_levels)) length(fit$batch_levels) else 1
  B_mod <- if (M > 1) B[-seq_len(M), , drop = FALSE] else B
  if (ncol(X_new) != nrow(B_mod))
    stop("X_new의 컬럼 수가 훈련 시 사용된 mod 디자인(B_mod)과 맞지 않습니다.")
  
  # --- 사전 계산 ---
  StS  <- crossprod(S)                                    # L x L
  Rinv <- chol2inv(chol(R + eps * diag(L)))               # L x L
  Dinv <- diag(1 / sqrt(pmax(sigma2_star, eps)), L, L)    # L x L
  SigmaAinv_star <- Dinv %*% Rinv %*% Dinv                # L x L
  
  # 공통 K, Q (배치불변 'star' 파라미터 사용)
  K <- SigmaAinv_star + StS / phi2_star + eps * diag(L)
  K <- (K + t(K)) / 2
  Q <- chol2inv(chol(K))                                  # L x L
  
  # --- A의 posterior mode (각 i에 대해 닫힌형) ---
  # A_i = [ X_i B_mod SigmaAinv_star + y_i S / phi2_star ] Q
  XB   <- X_new %*% B_mod                                 # n_new x L
  YS   <- (dat_new[, nonzero, drop = FALSE] %*% S) / phi2_star  # n_new x L
  A_new <- (XB %*% SigmaAinv_star + YS) %*% Q             # n_new x L
  
  # 선택적: 신호 재구성 (활성 엣지 차원만)
  Signal_new_nz <- A_new %*% t(S)                         # n_new x |nonzero|
  
  list(
    A = A_new,
    Signal_nz = Signal_new_nz,
    nonzero = nonzero,
    used = list(
      S = S,
      R = R,
      sigma2_star = sigma2_star,
      phi2_star = phi2_star,
      B_mod = B_mod,
      Q = Q,
      SigmaAinv_star = SigmaAinv_star
    )
  )
}