\name{SLACC}
\alias{SLACC}
\title{Sparse Latent Covariate-driven Connectome (SLACC) factorization}
\description{
Fits the SLACC model for symmetric brain connectivity data represented in vectorized upper-triangular form. The method factorizes each subject-specific connectivity matrix into a sum of sparse rank-1 symmetric patterns, with subject scores modeled as linear functions of covariates and batch/site membership. Estimation is carried out with a penalized expectation-maximization algorithm described in the manuscript, with optional sparsity regularization and annealing for updating latent loadings.
}
\usage{
SLACC(
  dat,
  mod = NULL,
  L = 5,
  batch = NULL,
  include_diag = TRUE,
  init = NULL,
  maxIter = 30,
  eps = 1e-3,
  U_maxIter = 200,
  U_eps = 1e-3,
  lambda_U = NULL,
  tau = NULL,
  annealing = TRUE,
  lambda_warmup_iter = NULL
)
}
\arguments{
  \item{dat}{An \eqn{n \times p} numeric matrix of vectorized connectivity data, where each row corresponds to one subject and each column corresponds to one upper-triangular entry of a symmetric \eqn{V \times V} connectivity matrix. The expected ordering is the same as used internally by the package's vectorization operator.}

  \item{mod}{An optional object specifying biological covariates to be preserved in the latent score model. In typical use this is a data frame or design-ready object passed to \code{model.matrix}. When \code{NULL}, the design includes only an intercept in the single-batch case and only batch indicators in the multi-batch case.}

  \item{L}{A positive integer giving the number of latent sparse rank-1 connectivity patterns. This is the latent dimension of the factorization. In the manuscript, \eqn{L} may be selected externally, for example by an information criterion.}

  \item{batch}{An optional vector or factor of length \code{n} indicating site or batch membership. If \code{NULL}, all subjects are treated as belonging to a single batch.}

  \item{include_diag}{Logical; if \code{TRUE}, diagonal entries are included in estimation. If \code{FALSE}, diagonal entries are excluded from the fitting objective and variance updates. For functional connectivity applications, diagonals are often omitted.}

  \item{init}{An optional list of initial values for the EM algorithm. If \code{NULL}, the function uses an internal higher-order SVD-based initialization. A custom initialization should be compatible with the internal parameterization and typically contains latent scores, loadings, regression coefficients, and variance components.}

  \item{maxIter}{Maximum number of EM iterations.}

  \item{eps}{Convergence tolerance for the outer EM algorithm, based on successive changes in the loading matrix after alignment.}

  \item{U_maxIter}{Maximum number of iterations for the inner loading-update routine. This applies to the ADMM or alternating least squares step used to update \code{U}.}

  \item{U_eps}{Convergence tolerance for the inner loading-update routine.}

  \item{lambda_U}{Nonnegative sparsity penalty level for the loading matrix \code{U}. If \code{NULL}, the default is \code{log(nrow(dat))}, following the manuscript's BIC-motivated penalty scale. Set to \code{0} to fit without sparsity penalization in the \code{U} update.}

  \item{tau}{Nonnegative truncation parameter used in the truncated lasso surrogate to the \eqn{L_0} penalty. If \code{NULL}, the default is \code{0.5 * sqrt(log(V * L) / n)}, matching the manuscript's default choice.}

  \item{annealing}{Logical; if \code{TRUE}, the sparsity penalty is gradually increased during early EM iterations before reaching \code{lambda_U}. This can improve numerical stability and reduce sensitivity to local minima.}

  \item{lambda_warmup_iter}{Optional integer giving the number of warm-up iterations used when \code{annealing = TRUE}. If \code{NULL}, an internal default is used.}
}
\details{
The fitted model represents each subject-specific connectivity matrix as a linear combination of sparse rank-1 symmetric patterns.
For subject \eqn{j} in batch/site \eqn{i}, the latent score for pattern \eqn{l} is modeled as a linear function of covariates plus a batch-specific random deviation,
which induces batch-specific latent variances and residual variances.

The function estimates the latent scores, sparse loading vectors, regression coefficients, and variance components jointly.
The implementation uses an EM algorithm, with a specialized update for the loading matrix \code{U}. When \code{lambda_U > 0},
a sparsity-inducing penalty is applied to the loading update; when \code{lambda_U = 0}, an unpenalized loading update is used.

The manuscript also discusses selection of \eqn{L} by an extended BIC. This function does not automatically search over multiple values of \code{L}; users should call \code{SLACC()} repeatedly if model comparison across candidate values of \code{L} is desired.
}
\value{
A list with three components:
\describe{
  \item{\code{estimates}}{A list of fitted model quantities, including latent subject scores \code{A}, vectorized latent patterns \code{S}, loading matrix \code{U}, regression coefficients \code{B}, latent covariance template \code{R}, batch-specific latent variances \code{sigma2}, batch-specific residual variances \code{phi2}, and residual means used by downstream harmonization routines.}
  \item{\code{input}}{A list recording model inputs and tuning parameters used during fitting, including the design matrix \code{X}, batch labels, \code{L}, penalty settings, convergence settings, and whether diagonals were included.}
  \item{\code{measure}}{A list of summary measures for the fitted model, including the log-likelihood, an effective degrees-of-freedom count, and the number of outer EM iterations completed.}
}
}
\references{
Zhang, R., Tuzhilina, E., & Park, J. Y. (2026). Sparse covariate-driven factorization of high-dimensional brain connectivity with application to site effect correction. arXiv preprint arXiv:2601.09525.
}
\examples{
set.seed(1)

n <- 20
V <- 5
p <- V * (V + 1) / 2
X <- data.frame(age = rnorm(n), sex = factor(sample(c("F", "M"), n, TRUE)))
batch <- factor(rep(c("site1", "site2"), each = n / 2))
dat <- matrix(rnorm(n * p), nrow = n, ncol = p)

fit <- SLACC(
  dat = dat,
  mod = X,
  L = 2,
  batch = batch,
  include_diag = TRUE,
  maxIter = 5,
  U_maxIter = 20
)

str(fit$estimates)
}
\keyword{models}
\keyword{multivariate}
