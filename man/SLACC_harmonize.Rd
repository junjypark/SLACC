\name{SLACC_harmonize}
\alias{SLACC_harmonize}
\title{Harmonize connectivity data using a fitted SLACC model}
\description{
Applies post-estimation harmonization based on a fitted \code{SLACC} object.
The routine removes batch-specific latent mean effects and rescales latent-factor and residual variability while preserving covariate-associated biological effects encoded in the fitted model.
}
\usage{
SLACC_harmonize(
  dat,
  fit,
  mod = NULL,
  batch = NULL
)
}
\arguments{
  \item{dat}{An \eqn{n \times p} numeric matrix of vectorized connectivity data to harmonize, with the same column convention used to fit the model. This may be the training data or new data from batches observed during fitting.}

  \item{fit}{A fitted object returned by \code{SLACC()}. The object must contain the estimated latent patterns, regression coefficients, and batch-specific variance parameters used for harmonization.}

  \item{mod}{Optional biological covariate information for \code{dat}. This should be compatible with the covariate specification used during model fitting. If omitted, the function uses the design structure implied by the fitted object.}

  \item{batch}{Optional batch/site labels for \code{dat}. These labels must correspond to batches seen in the training fit. If omitted, the function uses the training batch labels stored in \code{fit}.}
}
\details{
The manuscript defines harmonization as a two-part adjustment. First, batch-specific latent mean effects are removed while preserving the contribution of biological covariates. Second, batch-specific latent variances and residual variances are rescaled toward pooled targets.

In the manuscript, the pooled target latent variance for factor \eqn{l} is
\deqn{\sigma_l^{(h)} = \sqrt{\sum_i (n_i / n) \sigma_{il}^2}}
and the pooled residual scale is
\deqn{\phi^{(h)} = \sqrt{\sum_i (n_i / n) \phi_i^2}.}
Thus, the intended harmonization targets are weighted by batch sample sizes only.

For batch \eqn{i} and subject \eqn{j}, the adjusted latent scores take the form
\deqn{a_{ijl}^{(h)} = \frac{\sigma_l^{(h)}}{\hat\sigma_{il}}\left(\hat a_{ijl} - \hat\alpha_l - z_{ij}^\top \hat\theta_l - \hat\gamma_{il}\right) + \hat\alpha_l + z_{ij}^\top \hat\theta_l,}
with the residual term rescaled analogously by \eqn{\phi^{(h)} / \hat\phi_i}.
The harmonized connectivity vectors are then reconstructed from the adjusted latent scores and adjusted residuals.

To remain faithful to the manuscript, the biological covariates provided to \code{mod} should match those used when fitting \code{SLACC()}, and \code{batch} should identify only batches present in the training data.
}
\value{
A list containing harmonized representations derived from the fitted model. Depending on the package implementation, this typically includes harmonized connectivity data and may also include intermediate adjusted latent quantities used to reconstruct the harmonized data.
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

harm <- SLACC_harmonize(
  dat = dat,
  fit = fit,
  mod = X,
  batch = batch
)

str(harm)
}
\seealso{\code{\link{SLACC}}}
\keyword{models}
\keyword{multivariate}
