// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List bilinear_als3_cpp(const arma::cube& Y,     // V x V x n
                             const arma::cube& A,     // L x L x n
                             const arma::vec& w,      // n
                             arma::mat U,             // V x L
                             arma::mat V,             // V x L
                             const arma::mat& Qbar,   // L x L
                             const bool include_diag,
                             const int maxit,
                             const double tol,
                             const double gamma) {

  const int Vdim = Y.n_rows;
  const int Vdim2 = Y.n_cols;
  const int n = Y.n_slices;
  const int Ldim = A.n_rows;

  if (Vdim != Vdim2)
    stop("Y must be V x V x n");
  if ((int)A.n_rows != Ldim || (int)A.n_cols != Ldim || (int)A.n_slices != n)
    stop("A must be L x L x n, same n as Y");
  if ((int)w.n_elem != n)
    stop("w length must equal n");
  if ((int)U.n_rows != Vdim || (int)U.n_cols != Ldim)
    stop("U must be V x L");
  if ((int)V.n_rows != Vdim || (int)V.n_cols != Ldim)
    stop("V must be V x L");
  if ((int)Qbar.n_rows != Ldim || (int)Qbar.n_cols != Ldim)
    stop("Qbar must be L x L");

  double wsum = arma::accu(w);

  auto fro_sq = [](const arma::mat& M) {
    return accu(M % M);
  };

  for (int it = 0; it < maxit; ++it) {

    // ---------- U-step ----------
    arma::mat VtV = V.t() * V;        // L x L
    arma::mat GV  = arma::zeros<arma::mat>(Ldim, Ldim);
    arma::mat HV  = arma::zeros<arma::mat>(Vdim, Ldim);

    for (int i = 0; i < n; ++i) {
      arma::mat Ai = A.slice(i);     // L x L
      arma::mat Yi = Y.slice(i);     // V x V

      if (!include_diag) {
        arma::vec a = Ai.diag();                 // L
        arma::mat UV = U % V;                    // V x L
        arma::vec diag_fit = UV * a;             // V
        Yi.diag() = diag_fit;
      }

      GV += w(i) * (Ai * VtV * Ai.t());
      HV += w(i) * (Yi * V * Ai.t());
    }

    arma::mat K = GV + wsum * (VtV % Qbar) + gamma * arma::eye<arma::mat>(Ldim, Ldim);
    K = 0.5 * (K + K.t());

    arma::mat cholK;
    bool ok = arma::chol(cholK, K);
    if (!ok) {
      arma::mat K2 = K + 1e-8 * arma::eye<arma::mat>(Ldim, Ldim);
      ok = arma::chol(cholK, K2);
      if (!ok) stop("chol failed in U-step");
    }

    // U_new = solve(K, HV)
    arma::mat U_new = trans(arma::solve(trimatu(cholK),
                           arma::solve(trimatl(cholK.t()), HV.t())));

    // ---------- V-step ----------
    arma::mat UtU = U_new.t() * U_new;   // L x L
    arma::mat GU  = arma::zeros<arma::mat>(Ldim, Ldim);
    arma::mat HU  = arma::zeros<arma::mat>(Vdim, Ldim);

    for (int i = 0; i < n; ++i) {
      arma::mat Ai = A.slice(i);     // L x L
      arma::mat Yi = Y.slice(i);     // V x V

      if (!include_diag) {
        arma::vec a = Ai.diag();
        arma::mat UV = U_new % V;    // 여기서는 U_new, V
        arma::vec diag_fit = UV * a;
        Yi.diag() = diag_fit;
      }

      GU += w(i) * (Ai * UtU * Ai.t());
      HU += w(i) * (Yi.t() * U_new * Ai.t());
    }

    arma::mat K2 = GU + wsum * (UtU % Qbar) + gamma * arma::eye<arma::mat>(Ldim, Ldim);
    K2 = 0.5 * (K2 + K2.t());

    arma::mat cholK2;
    ok = arma::chol(cholK2, K2);
    if (!ok) {
      arma::mat K3 = K2 + 1e-8 * arma::eye<arma::mat>(Ldim, Ldim);
      ok = arma::chol(cholK2, K3);
      if (!ok) stop("chol failed in V-step");
    }

    arma::mat V_new = trans(arma::solve(trimatu(cholK2),
                           arma::solve(trimatl(cholK2.t()), HU.t())));

    double r_norm = std::sqrt(fro_sq(U_new - U) + fro_sq(V_new - V));

    double scale = std::sqrt( (double)Vdim * (double)Ldim );
    r_norm /= scale;

    U = U_new;
    V = V_new;

    if (r_norm < tol) break;
  }

  arma::mat U_final = 0.5 * (U + V);

  return Rcpp::List::create(
    Rcpp::Named("U") = U_final,
    Rcpp::Named("V") = U_final
  );
}