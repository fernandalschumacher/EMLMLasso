#' Fit a lasso-penalized linear mixed model via EM, for a single lambda
#'
#' Internal engine, called by [EMLMLasso()] both for a single-lambda fit and
#' for every grid point of a BIC search. Not exported -- use [EMLMLasso()].
#'
#' Fits a linear mixed model with a lasso penalty on the fixed-effects
#' coefficients for a single, fixed value of the tuning parameter `lambda`,
#' using an Expectation-Maximization (EM) algorithm. At each E-step, the random
#' effects are treated as missing data and their conditional means and
#' covariances are computed in closed form; at each M-step, the fixed
#' effects are updated by a coordinate-descent lasso fit ([glmnet::glmnet()])
#' on the random-effects-adjusted response, and the variance components are
#' updated in closed form.
#'
#' Following Oliveira et al. (2026), the model is
#' \eqn{Y_i = X_i \beta + Z_i b_i + \epsilon_i}, with
#' \eqn{b_i \sim N_q(0, D)} independent of \eqn{\epsilon_i \sim N_{n_i}(0,
#' \sigma^2 I_{n_i})}, so the covariance for subject/cluster
#' \eqn{i} is \eqn{\Sigma_i = Z_i D Z_i' + \sigma^2 I_{n_i}}.
#'
#' `lambda` here plays the role of \eqn{\lambda^{(k)}} in the M-step of
#' Algorithm 1 of Oliveira et al. (2026) (their Equation 8) -- the penalty
#' passed directly to [glmnet::glmnet()]. Unlike their Algorithm 1, `lambda`
#' is fixed for the duration of one call rather than re-selected by
#' cross-validation at every EM iteration; see [EMLMLasso()] and `NEWS.md`
#' for how this package selects it instead. Following Equation 1 of the same
#' paper, no intercept is fit by default (see `intercept`).
#'
#' @param data a data.frame containing all variables referenced by
#'   `formFixed`, `formRandom` and `groupVar`.
#' @param formFixed two-sided formula for the fixed effects, already resolved
#'   (no `y ~ .` shorthand -- that is expanded by [EMLMLasso()] before this is
#'   called).
#' @param groupVar character string naming the subject/cluster variable in `data`.
#' @param formRandom one-sided formula for the random effects.
#' @param lambda a single non-negative penalty value.
#' @param intercept if `TRUE`, force an intercept into the fixed-effects
#'   design matrix; if `FALSE` (the default), force it out.
#' @param center if `TRUE` and `intercept = FALSE`, centre the response and
#'   every fixed-effect covariate at their overall means before fitting.
#' @param precision convergence tolerance on the Euclidean distance between
#'   successive parameter vectors (beta, sigma2, D).
#' @param MaxIter maximum number of EM iterations.
#' @param penalty_factor optional vector of per-coefficient penalty
#'   multipliers, in the same order as the columns of the fixed-effects
#'   design matrix. If `NULL`, all coefficients are penalized equally except
#'   an intercept, if present, which is left unpenalized.
#' @param standardize_x if `TRUE`, each fixed-effect covariate (other than an
#'   intercept, if present) is scaled to unit variance once before fitting,
#'   and coefficients are converted back to the original scale of `x` before
#'   being returned.
#' @param refit if `TRUE`, after the lasso fit converges, also fits an
#'   ordinary (unpenalized) `lme4::lmer()` model on just the covariates the
#'   lasso left non-zero (same random-effects structure), to correct the
#'   shrinkage bias lasso point estimates always carry. Requires the `lme4`
#'   package.
#' @param verbose if `TRUE`, prints one progress line every 10 iterations.
#'
#' @return a plain list with elements `coefficients`, `sigma2`, `D`, `ranef`,
#'   `lambda`, `penalty_factor`, `loglik`, `AIC`, `BIC`, `npar`, `iter`,
#'   `converged`, `glmnet_converged`, `n`, `N`, `x`, `z`, `y`, `nj`, `groups`,
#'   `call`, `formFixed`, `formRandom`, `groupVar`, and, if `refit = TRUE`,
#'   `refit_fit`, `refit_coefficients`, `refit_formula`, `refit_singular` --
#'   see [EMLMLasso()]'s `@return` for what each of these means; this
#'   internal engine's output becomes exactly that user-facing object's
#'   fields when a single `lambda` is requested, or `best_model`'s fields
#'   when it is called from a BIC search.
#' @noRd
EMLasso_LMM <- function(data, formFixed, groupVar, formRandom = ~1, lambda,
                         intercept = FALSE, center = TRUE, precision = 1e-6, MaxIter = 200,
                         penalty_factor = NULL, standardize_x = TRUE,
                         refit = FALSE, verbose = interactive()) {
  cl <- match.call()
  if (missing(lambda) || length(lambda) != 1L || !is.numeric(lambda) || lambda < 0) {
    stop("`lambda` must be a single non-negative number.", call. = FALSE)
  }
  if (isTRUE(refit) && !requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for `refit = TRUE`. Install it with ",
         "install.packages(\"lme4\").", call. = FALSE)
  }

  dp <- prepare_lmm_data(formFixed, formRandom, groupVar, data, intercept = intercept)
  y <- dp$y; x <- dp$x; z <- dp$z; nj <- dp$nj

  n <- length(nj); N <- sum(nj); p <- ncol(x); q <- ncol(z)

  if (is.null(penalty_factor)) {
    penalty_factor <- default_penalty_factor(x)
  } else if (length(penalty_factor) != p) {
    stop("`penalty_factor` must have length ", p, " (one entry per fixed-effect ",
         "coefficient, including the intercept if present).", call. = FALSE)
  }

  has_intercept <- "(Intercept)" %in% colnames(x)
  x_orig <- x; y_orig <- y

  # --- Centre y and x (only without an intercept column, which would
  # otherwise absorb the means) ---
  centered <- isTRUE(center) && !has_intercept
  if (centered) {
    y_mean <- mean(y)
    x_means <- colMeans(x)
    y <- y - y_mean
    x <- sweep(x, 2, x_means, "-")
  }

  # --- Standardize x (scale only; centring, if any, is done above) ---
  if (standardize_x) {
    x_sd <- apply(x, 2, stats::sd)
    if (has_intercept) x_sd[["(Intercept)"]] <- 1
    x_sd[x_sd == 0] <- 1 # guard against a constant/degenerate covariate
    x <- sweep(x, 2, x_sd, "/")
  }
  glmnet_standardize <- !standardize_x

  # --- Initialization ---
  betaI <- glmnet::glmnet(x, y, lambda = lambda, family = "gaussian",
                          intercept = FALSE, alpha = 1, penalty.factor = penalty_factor,
                          standardize = glmnet_standardize)
  beta <- betaI$beta[, 1]

  sigmae <- sum((y - x %*% beta)^2) / N
  sigmae <- max(sigmae, 1e-4) # variance floor

  D <- diag(q)
  loglik <- lmm_loglik(beta, sigmae, D, y, x, z, nj)
  loglikvec <- numeric(MaxIter)
  loglikvec[1] <- loglik
  teta <- c(beta, sigmae, D[upper.tri(D, diag = TRUE)])

  criterio <- 1; count <- 0
  end_indices <- cumsum(nj)
  start_indices <- c(1, end_indices[-n] + 1)
  glmnet_converged <- TRUE

  while (criterio > precision && count < MaxIter) {
    count <- count + 1
    sum3 <- 0; sum4 <- matrix(0, q, q)
    b <- matrix(0, q, n); zm <- matrix(0, nrow = N, ncol = 1)

    # --- E-step ---
    for (j in seq_len(n)) {
      idx <- start_indices[j]:end_indices[j]
      y1 <- y[idx]
      x1 <- x[idx, , drop = FALSE]
      z1 <- z[idx, , drop = FALSE]

      med <- x1 %*% beta
      Psi <- z1 %*% D %*% t(z1) + sigmae * diag(nj[j])
      Lambda <- solve(solve(D) + t(z1) %*% z1 / sigmae, tol = 1e-30)

      b[, j] <- Lambda %*% t(z1) %*% (y1 - med) / sigmae
      y1m <- y1 - z1 %*% b[, j]

      sum3 <- sum3 + t(y1m - med) %*% (y1m - med) + sum(diag(Lambda %*% t(z1) %*% z1))
      sum4 <- sum4 + Lambda + b[, j] %*% t(b[, j])
      zm[idx, 1] <- y1m
    }

    # --- M-step ---
    la.eq <- do.call(glmnet::glmnet, c(
      list(x = x, y = zm, lambda = lambda, family = "gaussian",
           intercept = FALSE, alpha = 1, penalty.factor = penalty_factor,
           standardize = glmnet_standardize),
      glmnet_thresh_arg(1e-8)
    ))

    if (isTRUE(la.eq$jerr != 0)) {
      warning("glmnet's coordinate-descent lasso did not converge in the ",
              "M-step at iteration ", count, " (lambda = ", signif(lambda, 4),
              "); returning the fit from the last iteration that did ",
              "converge instead of a corrupted one.", call. = FALSE)
      glmnet_converged <- FALSE
      count <- count - 1
      break
    }

    beta <- la.eq$beta[, 1]

    sigmae <- max(as.numeric(sum3) / N, 1e-4)
    D <- sum4 / n

    if (verbose && count %% 10 == 0) {
      cat("  Iter:", count, "| Sigma_e2:", round(sigmae, 6),
          "| Lambda:", round(lambda, 6), "\n")
    }

    # --- Convergence ---
    loglik <- lmm_loglik(beta, sigmae, D, y, x, z, nj)
    if (is.infinite(loglik) || is.nan(loglik)) break

    loglikvec[count] <- loglik
    teta1 <- teta
    teta <- c(beta, sigmae, D[upper.tri(D, diag = TRUE)])
    criterio <- sqrt(sum((teta1 - teta)^2))

  }

  teta <- c(beta, sigmae, D[upper.tri(D, diag = TRUE)])
  npar <- length(teta) - sum(beta == 0)
  AIC <- -2 * loglik + 2 * npar
  BIC <- -2 * loglik + log(N) * npar

  # --- Undo the internal standardization and centring ---
  if (standardize_x) beta <- beta / x_sd
  x <- x_orig; y <- y_orig
  implied_intercept <- if (centered) y_mean - sum(x_means * beta) else 0

  ranef <- t(b)
  rownames(ranef) <- dp$groups
  colnames(ranef) <- colnames(z)
  names(beta) <- colnames(x)

  obj <- list(
    coefficients = beta,
    implied_intercept = implied_intercept,
    centered = centered,
    sigma2 = sigmae,
    D = D,
    ranef = ranef,
    lambda = lambda,
    penalty_factor = penalty_factor,
    loglik = loglikvec[seq_len(count)],
    AIC = AIC, BIC = BIC, npar = npar,
    iter = count, converged = (criterio <= precision) && glmnet_converged,
    glmnet_converged = glmnet_converged,
    n = n, N = N,
    x = x, z = z, y = y, nj = nj, groups = dp$groups,
    call = cl, formFixed = dp$formFixed,
    formRandom = dp$formRandom, groupVar = groupVar
  )

  if (isTRUE(refit)) {
    refit_res <- refit_lmm_impl(obj)
    obj$refit_fit <- refit_res$refit_fit
    obj$refit_coefficients <- refit_res$refit_coefficients
    obj$refit_formula <- refit_res$refit_formula
    obj$refit_singular <- refit_res$refit_singular
  }

  obj
}
