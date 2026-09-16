#' Simulate longitudinal data from a sparse linear mixed model
#'
#' Generates balanced longitudinal/clustered data from a linear mixed model
#' with `p` fixed-effect covariates (only the first `pstar` have a non-zero
#' coefficient) and a random intercept + random slope on `time`. This is the
#' data-generating mechanism used in the package examples, tests and
#' vignette (adapted from the original simulation study for this method).
#'
#' @param n number of subjects/clusters.
#' @param nj number of observations per subject (balanced design).
#' @param p total number of fixed-effect covariates.
#' @param pstar number of covariates with a non-zero true coefficient
#'   (`beta_value`); the remaining `p - pstar` have a true coefficient of 0.
#'   Must be `<= p`.
#' @param beta_value the common value used for the non-zero true coefficients.
#' @param sigmae true residual variance.
#' @param D true 2x2 covariance matrix of the random intercept and random
#'   slope on `time`. Defaults to `matrix(c(1, 0.25, 0.25, 1), 2, 2)`.
#' @param seed optional seed passed to [set.seed()] for reproducibility.
#'
#' @return A data.frame with columns `id`, `time`, `y`, and `x1, ..., xp`.
#'
#' @examples
#' dat <- simulate_lasso_lmm(n = 20, nj = 4, p = 6, pstar = 2, seed = 42)
#' head(dat)
#'
#' @export
simulate_lasso_lmm <- function(n = 50, nj = 5, p = 10, pstar = 4, beta_value = 1,
                                sigmae = 1, D = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (pstar > p) stop("`pstar` cannot exceed `p`.", call. = FALSE)
  if (is.null(D)) D <- matrix(c(1, 0.25, 0.25, 1), 2, 2)

  nj_vec <- rep(nj, n)
  N <- sum(nj_vec)
  end_indices <- cumsum(nj_vec)
  start_indices <- c(1, end_indices[-n] + 1)

  X <- matrix(stats::rnorm(N * p), ncol = p)
  beta <- c(rep(beta_value, pstar), rep(0, p - pstar))
  time <- unlist(lapply(nj_vec, seq_len))
  id <- rep(seq_len(n), times = nj_vec)
  z <- cbind(1, time)

  y <- numeric(N)
  for (j in seq_len(n)) {
    idx <- start_indices[j]:end_indices[j]
    x1 <- X[idx, , drop = FALSE]
    z1 <- z[idx, , drop = FALSE]
    bj <- mvtnorm::rmvnorm(1, rep(0, 2), D)
    ej <- mvtnorm::rmvnorm(1, rep(0, nj_vec[j]), sigmae * diag(nj_vec[j]))
    y[idx] <- x1 %*% beta + z1 %*% t(bj) + t(ej)
  }

  colnames(X) <- paste0("x", seq_len(p))
  data.frame(id = id, time = time, y = as.numeric(y), X, check.names = FALSE)
}
