# SLACC

**SLACC** implements sparse latent covariate-driven factorization for high-dimensional symmetric brain connectivity data, with an application to site/batch effect correction in multi-site studies. It is designed for brain connectivity matrices such as functional connectivity (FC). The underlying model factorizes each subject-specific connectivity matrix into a sum of sparse rank-1 symmetric patterns, with subject-level latent scores modeled by biological covariates and batch/site effects.

## Summary

For subject `j` from site `i`, the manuscript models a symmetric connectivity matrix `Y_ij` as

`Y_ij = sum_{l=1}^L a_ijl u_l u_l^T + E_ij`,

where:

- `u_l` is a sparse loading vector defining the `l`th latent connectivity pattern,
- `a_ijl` is the subject-specific score for pattern `l`,
- `E_ij` is a site-specific residual term.

The latent scores are modeled as

`a_ijl = x_ij^T beta_l + delta_ijl`,

allowing biological covariates and site effects to enter through the latent representation. In the manuscript, site-specific heterogeneity is represented through:

- site-specific latent mean effects,
- site-specific latent variances `sigma^2_il`, and
- site-specific residual variances `phi^2_i`.

Model fitting is carried out using a penalized EM algorithm. Sparsity in the loading matrix is induced through a truncated-lasso surrogate to an `L0` penalty, and the implementation supports annealing during the loading update to improve numerical stability.

## Main functions

### `SLACC()`

Fits the sparse latent covariate-driven factorization model.

Typical inputs are:

- `dat`: an `n x p` numeric matrix containing vectorized upper-triangular entries of symmetric connectivity matrices,
- `mod`: optional biological covariates,
- `L`: number of latent sparse rank-1 connectivity patterns,
- `batch`: optional site/batch labels.

See `?SLACC` for the full argument list.

### `SLACC_harmonize()`

Applies post-estimation harmonization using a fitted `SLACC` object.

The intended harmonization described in the manuscript:

1. removes batch-specific latent mean effects while preserving biological covariate effects,
2. rescales batch-specific latent variances toward pooled latent targets,
3. rescales batch-specific residual variances toward a pooled residual target,
4. reconstructs harmonized connectivity data from the adjusted latent scores and residuals.

See `?SLACC_harmonize` for details.

## Installation

The package can be installed from GitHub with `remotes`:

```r
install.packages("remotes")
remotes::install_github("junjypark/SLACC")
```

## Basic workflow

```r
library(SLACC)

## dat: n x p matrix of vectorized symmetric connectivities
## mod: data frame or design-ready covariates to preserve
## batch: site/scanner labels

fit <- SLACC(
  dat = dat,
  mod = mod,
  L = 5,
  batch = batch,
  include_diag = TRUE,
  maxIter = 30,
  eps = 1e-3,
  U_maxIter = 200,
  U_eps = 1e-3,
  lambda_U = NULL,
  tau = NULL,
  annealing = TRUE
)

harm <- SLACC_harmonize(
  dat = dat,
  fit = fit,
  mod = mod,
  batch = batch
)
```

## Example with simulated data

The codes below mirror the basic simulation workflow used in the manuscript: specify sparse latent patterns `U`, covariates `X`, site-specific latent variances, site-specific residual variances, generate `dat`, and fit `SLACC()`. 

```r
library(SLACC)
library(MASS)

set.seed(1)

## Example objects typically prepared in advance
U <- readRDS("U_setting2.rds")
X <- readRDS("X_n100.rds")
beta <- readRDS("beta.rds")

n <- nrow(X)
L <- ncol(U)
V <- nrow(U)
p <- V * (V + 1) / 2

sigma2_g1 <- 8:12
sigma2_g2 <- 12:8
phi2 <- c(0.8, 1.2)
batch <- c(rep("a", n / 2), rep("b", n / 2))

## Ltrans() should map a symmetric V x V matrix to its upper-triangular vectorization
M <- 0
for (l in 1:L) {
  S_l <- Ltrans(tcrossprod(U[, l]))
  A_l <- c(
    -0.3 + X[1:(n / 2), ] %*% beta[, l] + rnorm(n / 2, 0, sqrt(sigma2_g1[l])),
     0.3 + X[-(1:(n / 2)), ] %*% beta[, l] + rnorm(n / 2, 0, sqrt(sigma2_g2[l]))
  )
  M <- M + tcrossprod(A_l, S_l)
}

E <- matrix(NA_real_, n, p)
E[1:(n / 2), ] <- matrix(rnorm((n / 2) * p, 0, sqrt(phi2[1])), n / 2, p)
E[-(1:(n / 2)), ] <- matrix(rnorm((n / 2) * p, 0, sqrt(phi2[2])), n / 2, p)

dat <- M + E

fit <- SLACC(
  dat = dat,
  mod = X,
  L = L,
  batch = batch
)
```

If the goal is post-hoc harmonization after fitting, use:

```r
harm <- SLACC_harmonize(
  dat = dat,
  fit = fit,
  mod = X,
  batch = batch
)
```

## Input format

`SLACC` expects a matrix of vectorized symmetric connectivity data.

- If the original connectivity matrix is `V x V`, then the input has `p = V(V + 1)/2` columns when diagonal entries are included.
- Each row corresponds to one subject.
- Columns should follow a consistent upper-triangular vectorization convention.

For functional connectivity applications, users may choose whether to include diagonal entries in fitting via `include_diag`.

## Output

A fitted `SLACC` object typically contains:

- `estimates`: fitted latent scores, loading vectors, regression coefficients, and variance parameters,
- `input`: model inputs and tuning settings,
- `measure`: fit summaries such as log-likelihood and iteration counts.

The harmonization function returns harmonized representations derived from the fitted model and may also return intermediate latent quantities used during reconstruction.

## Notes

- The manuscript discusses selecting the latent dimension `L` using an extended BIC across candidate values of `L`. The exported fitting function takes `L` as an input; users can compare multiple fits externally.
- Biological covariates supplied to `mod` should represent effects to preserve during harmonization.
- Batch labels passed to `SLACC_harmonize()` should correspond to batches represented in the fitted model.
- The simulation-style example above is intentionally simplified and does not include parallelization or repeated runs.

## Documentation

Package help pages:

- `?SLACC`
- `?SLACC_harmonize`

## Reference

Zhang, R., Tuzhilina, E., & Park, J. Y. (2026). Sparse covariate-driven factorization of high-dimensional brain connectivity with application to site effect correction. arXiv preprint arXiv:2601.09525.
