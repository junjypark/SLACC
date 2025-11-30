prepare_elements = function(dat, A, U, L, phi2, tau=NULL, ni){
  phi2 = rep(unlist(phi2), ni)
  n = nrow(dat)
  p = ncol(dat)
  V = (sqrt(1+8*p)-1)/2
  Y = array(NA, dim=c(n, V, V))
  X = array(NA, dim=c(n, L, L))
  for (i in 1:n){
    Y[i,,] = Ltrinv(dat[i,],V)
    X[i,,] = diag(A[i,], ncol=L)
  }
  subj_wts = 1/phi2
  return(list( Y=Y, X=X, subj_wts=subj_wts, B_wts=B_wts ))
}