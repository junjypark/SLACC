// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// sign 
inline double sign_double(double x) {
  if (x > 0) return 1.0;
  if (x < 0) return -1.0;
  return 0.0;
}

// [[Rcpp::export]]
Rcpp::List bilinear_admm_cpp(const arma::cube &Ycube,   // V x V x n
                              const arma::cube &Acube,   // L x L x n
                              const arma::vec  &w,       // n
                              const arma::mat  &C,       // V x L
                              const double      lambda, 
                              arma::mat         U,       // V x L
                              arma::mat         V,       // V x L
                              const arma::mat  &Qbar,    // L x L
                              const bool        include_diag = true,
                              const double      rho   = 1.0,
                              const double      eta   = 1.0,
                              const double      gamma_par = 1e-6,
                              const int         maxit = 100,
                              const double      tol   = 1e-3,
                              const int         snap_to_Z_every = 1) {
  
  int Vdim = Ycube.n_rows;
  int n = Ycube.n_slices;
  int Ldim = U.n_cols;
  
  if (Ycube.n_cols != Vdim || (int)w.n_elem != n) {
    stop("Dimension mismatch in Ycube / w");
  }
  if (Acube.n_rows != (unsigned)Ldim || Acube.n_cols != (unsigned)Ldim || (int)Acube.n_slices != n) {
    stop("Dimension mismatch in Acube");
  }
  if ( (int)C.n_rows != Vdim || (int)C.n_cols != Ldim ) {
    stop("C must be V x L");
  }
  if ( (int)Qbar.n_rows != Ldim || (int)Qbar.n_cols != Ldim ) {
    stop("Qbar must be L x L");
  }
  
  double wsum = arma::accu(w);
  
  arma::mat ZU = U;
  arma::mat ZV = V;
  arma::mat WU(Vdim, Ldim, fill::zeros);
  arma::mat WV(Vdim, Ldim, fill::zeros);
  arma::mat Yeq(Vdim, Ldim, fill::zeros);
  
  auto fro_sq = [](const arma::mat &M) {
    return arma::accu(arma::square(M));
  };
  
  arma::mat UV(Vdim, Ldim);
  
  for (int it = 0; it < maxit; ++it) {
    
    // ---------- U-step ----------
    arma::mat VtV = V.t() * V;           // L x L
    arma::mat GV(Ldim, Ldim, fill::zeros);
    arma::mat HV(Vdim, Ldim, fill::zeros);
    
    for (int i = 0; i < n; ++i) {
      arma::mat Ai = Acube.slice(i);    // L x L
      arma::mat Yi = Ycube.slice(i);    // V x V
      
      if (!include_diag) {
        // diag(Yi)를 현재 (U, V, a_ij)로부터의 fitted diag로 대체
        arma::vec a = Ai.diag();    // L
        UV = U % V;                 // V x L (elementwise)
        arma::vec diag_fit = UV * a;   // V x 1
        Yi.diag() = diag_fit;
      }
      
      GV += w(i) * (Ai * VtV * Ai.t());        // L x L
      HV += w(i) * (Yi * V * Ai.t());          // V x L
    }
    
    arma::mat K = GV 
      + (rho + eta + gamma_par) * arma::eye<arma::mat>(Ldim, Ldim)
      + wsum * (VtV % Qbar);                  // Hadamard with Qbar
    
    K = 0.5 * (K + K.t());
    K += 1e-8 * arma::eye<arma::mat>(Ldim, Ldim);
    
    arma::mat T_U = rho * (ZU - WU) + eta * (V - Yeq);
    arma::mat RHS = HV + T_U;                 // V x L
    
    // solve K * U^T = RHS^T
    arma::mat Ut = arma::solve(K, RHS.t(), arma::solve_opts::fast);
    U = Ut.t();                               // V x L
    
    // ZU soft-threshold + dual
    arma::mat ZU_old = ZU;
    arma::mat tmpU   = U + WU;
    arma::mat shrinkU = arma::abs(tmpU) - (lambda / rho) * C;
    shrinkU.transform( [](double x){ return (x > 0.0) ? x : 0.0; } );
    ZU = arma::sign(tmpU) % shrinkU;
    WU += (U - ZU);
    
    // ---------- V-step ----------
    arma::mat UtU = U.t() * U;            // L x L
    arma::mat GU(Ldim, Ldim, fill::zeros);
    arma::mat HU(Vdim, Ldim, fill::zeros);
    
    for (int i = 0; i < n; ++i) {
      arma::mat Ai = Acube.slice(i);     // L x L
      arma::mat Yi = Ycube.slice(i);     // V x V
      
      if (!include_diag) {
        arma::vec a = Ai.diag();
        UV = U % V;
        arma::vec diag_fit = UV * a;
        Yi.diag() = diag_fit;
      }
      
      GU += w(i) * (Ai * UtU * Ai.t());
      HU += w(i) * (Yi.t() * U * Ai.t());
    }
    
    K = GU 
      + (rho + eta + gamma_par) * arma::eye<arma::mat>(Ldim, Ldim)
      + wsum * (UtU % Qbar);
    
    K = 0.5 * (K + K.t());
    K += 1e-8 * arma::eye<arma::mat>(Ldim, Ldim);
    
    arma::mat T_V = rho * (ZV - WV) + eta * (U + Yeq);
    RHS = HU + T_V;                          // V x L
    
    arma::mat Vt = arma::solve(K, RHS.t(), arma::solve_opts::fast);
    V = Vt.t();                              // V x L
    
    // ZV soft-threshold + dual
    arma::mat ZV_old = ZV;
    arma::mat tmpV   = V + WV;
    arma::mat shrinkV = arma::abs(tmpV) - (lambda / rho) * C;
    shrinkV.transform( [](double x){ return (x > 0.0) ? x : 0.0; } );
    ZV = arma::sign(tmpV) % shrinkV;
    WV += (V - ZV);
    
    // equality dual (U = V)
    Yeq += (U - V);
    
    // snap
    if (snap_to_Z_every > 0 && ((it + 1) % snap_to_Z_every == 0)) {
      U = ZU;
      V = ZV;
    }
    
    double r_norm = std::sqrt(
      fro_sq(U - ZU) + fro_sq(V - ZV) + fro_sq(U - V)
    );
    double s_norm = rho * std::sqrt(
      fro_sq(ZU - ZU_old) + fro_sq(ZV - ZV_old)
    );
    
    if (std::max(r_norm, s_norm) < tol) {
      // Rcpp::Rcout << "ADMM converged at iter " << (it+1) << std::endl;
      break;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("U") = ZU,
    Rcpp::Named("V") = ZV,
    Rcpp::Named("rho") = rho
  );
}