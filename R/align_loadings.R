align_loadings <- function(U, Uhat, method = c("cosine", "corr")) {
  method <- match.arg(method)
  
  stopifnot(nrow(U) == nrow(Uhat))
  k  <- ncol(U)
  kh <- ncol(Uhat)
  if (kh < k) stop("Uhat has fewer columns than U; increase rank of Uhat or reduce k.")
  
  # column-normalize (so t(U) %*% Uhat is cosine similarity)
  norm_cols <- function(M) {
    s <- sqrt(colSums(M^2)); s[s == 0] <- 1
    sweep(M, 2, s, "/")
  }
  U_n   <- norm_cols(U)
  Uh_n  <- norm_cols(Uhat)
  
  # similarity matrix S (k x kh)
  if (method == "cosine") {
    S <- t(U_n) %*% Uh_n
  } else {
    # Pearson correlation per column pair (same as cosine after centering,
    # but included for completeness)
    center_cols <- function(M) sweep(M, 1, rowMeans(M), "-")
    U_c  <- center_cols(U); Uh_c <- center_cols(Uhat)
    U_c  <- norm_cols(U_c); Uh_c <- norm_cols(Uh_c)
    S <- t(U_c) %*% Uh_c
  }
  
  # Solve assignment on cost = 1 - |S| (maximize |S|)
  use_clue <- requireNamespace("clue", quietly = TRUE)
  if (use_clue) {
    # Pad to square if kh > k (common when Uhat has extra components)
    if (kh > k) {
      pad <- matrix(0, nrow = kh - k, ncol = kh) # zeros so they won’t be chosen
      cost <- rbind(1 - abs(S), pad)
    } else {
      cost <- 1 - abs(S)
    }
    assignment <- clue::solve_LSAP(cost) # returns length = nrow(cost)
    perm_full  <- as.integer(assignment)
    perm <- perm_full[seq_len(k)]        # first k match U’s columns
  } else {
    # Greedy fallback (not guaranteed optimal, but works without extra packages)
    warning("Package 'clue' not found; using greedy matching. Install 'clue' for optimal assignment.")
    perm <- integer(k)
    used <- rep(FALSE, kh)
    for (i in seq_len(k)) {
      j <- which.max((1 - used) * abs(S[i, ]))
      perm[i] <- j
      used[j] <- TRUE
    }
  }
  
  # Determine signs so that matched similarity is nonnegative
  signed_sim <- mapply(function(i, j) S[i, j], seq_len(k), perm)
  signs <- ifelse(signed_sim >= 0, 1, -1)
  
  # Build aligned Uhat: Uhat_aligned[, i] = signs[i] * Uhat[, perm[i]]
  Uhat_aligned <- Uhat[, perm, drop = FALSE] %*% diag(signs, nrow = k, ncol = k)
  
  list(
    Uhat_aligned = Uhat_aligned,  # columns ordered like U, signs fixed
    permutation  = perm,          # for each i (column of U), matched column index in Uhat
    signs        = signs,         # +1 / -1 applied to matched columns
    similarity   = S,             # raw similarity matrix (before abs/sign)
    matched_sim  = signed_sim,    # signed similarities after matching
    alignment_score = sum(abs(signed_sim)) / k  # average |cosine| across matched columns
  )
}