#' Fit a lasso-penalized linear mixed model via EM
#'
#' Fits a linear mixed model with a lasso penalty on the fixed-effects
#' coefficients for a single, fixed value of the tuning parameter `lambda`,
#' using an Expectation-Maximization (EM) algorithm. At each E-step, random
#' effects are integrated out in closed form; at each M-step, the fixed
#' effects are updated by a coordinate-descent lasso fit ([glmnet::glmnet()])
#' on the random-effects-adjusted response, and the variance components are
#' updated in closed form.
#'
#' Following Oliveira et al. (2026), the model is
#' \eqn{Y_i = X_i \beta + Z_i b_i + \epsilon_i}, with
#' \eqn{b_i \sim N_q(0, D)} independent of \eqn{\epsilon_i \sim N_{n_i}(0,
#' \sigma^2 I_{n_i})}, so the covariance for subject/cluster
#' \eqn{i} is \eqn{\Psi_i = Z_i D Z_i' + \sigma^2 I_{n_i}}.
#'
#' `lambda` here plays the role of \eqn{\lambda^{(k)}} in the M-step of
#' Algorithm 1 of Oliveira et al. (2026) (their Equation 8) -- the penalty
#' passed directly to [glmnet::glmnet()]. Unlike their Algorithm 1, `lambda`
#' is now fixed for the duration of one call rather than re-selected by
#' cross-validation at every EM iteration; see [EMLasso_BIC_Search()] and
#' `NEWS.md` for how this package selects it instead. Following Equation 1
#' of the same paper, no intercept is fit by default (see `intercept`).
#'
#' @references
#' Oliveira DCR, Schumacher FL, Lachos VH (2026). The use of the EM algorithm
#' for regularization problems in high-dimensional linear mixed-effects
#' models. *Statistical Methods in Medical Research*, 35(2), 240-254.
#' \doi{10.1177/09622802251399913}
#'
#' @param data a data.frame containing all variables referenced by
#'   `formFixed`, `formRandom` and `groupVar`.
#' @param formFixed two-sided formula for the fixed effects, e.g. `y ~ x1 + x2`.
#'   `y ~ .` is also accepted, meaning every other column of `data` is a
#'   fixed-effect candidate except the response and `groupVar` itself.
#' @param groupVar character string naming the subject/cluster variable in `data`.
#' @param formRandom one-sided formula for the random effects, e.g. `~1`
#'   (random intercept, the default) or `~time` (random intercept and slope).
#' @param lambda a single non-negative penalty value.
#' @param intercept if `TRUE`, force an intercept into the fixed-effects
#'   design matrix; if `FALSE` (the default), force it out. This always wins
#'   over whatever `formFixed` itself implies -- `y ~ x1 + x2` and
#'   `y ~ x1 + x2 - 1` are treated identically.
#' @param precision convergence tolerance on the Euclidean distance between
#'   successive parameter vectors (beta, sigma2, D). Default `1e-6`.
#' @param MaxIter maximum number of EM iterations. Default `200`.
#' @param penalty_factor optional vector of per-coefficient penalty
#'   multipliers, in the same order as the columns of the fixed-effects
#'   design matrix (see [stats::model.matrix()] applied to `formFixed`). If
#'   `NULL` (the default), all coefficients are penalized equally except an
#'   intercept, if present, which is left unpenalized. Set an entry to `0`
#'   to force a coefficient to stay in the model.
#' @param standardize_x if `TRUE` (the default), each fixed-effect covariate
#'   (other than an intercept, if `formFixed` includes one) is scaled to
#'   unit variance once before fitting, and coefficients are converted back
#'   to the original scale of `x` before being returned -- `coef()`,
#'   `fitted()`, etc. are unaffected by this option. This is
#'   similar to what `glmnet::glmnet()`'s own default
#'   (`standardize = TRUE`) already does internally at every call, but
#'   doing it once, up front, avoids recomputing it at every M-step of every
#'   EM iteration, and gives noticeably better numerical conditioning for
#'   many correlated, high-dimensional covariates. Set to `FALSE` to fit on
#'   the original scale of `x` instead (relying only on glmnet's own default
#'   internal standardization, recomputed every M-step).
#' @param refit if `TRUE`, after the lasso fit converges, also fits an
#'   ordinary (unpenalized) `lme4::lmer()` model on just the covariates the
#'   lasso left non-zero (same random-effects structure), to correct the
#'   shrinkage bias lasso point estimates always carry. Requires the `lme4`
#'   package. **The refit's own standard errors, p-values and confidence
#'   intervals must not be used or reported** as they treat the selected
#'   covariate set as fixed in advance, when it was in fact chosen using
#'   this same data, so they do not account for that selection step and
#'   understate the true uncertainty (classic post-selection inference).
#'   The original lasso fit (`coefficients`, `AIC`, `BIC`,
#'   etc.) is always kept unchanged alongside it. Default `FALSE`.
#' @param verbose if `TRUE`, prints one progress line every 10 iterations.
#'   Defaults to `interactive()`.
#'
#' @return An object of class `"EMLasso_LMM"`, which is a list with elements:
#'   \describe{
#'     \item{coefficients}{named vector of estimated fixed-effects coefficients.}
#'     \item{sigma2}{estimated residual variance.}
#'     \item{D}{estimated random-effects covariance matrix.}
#'     \item{ranef}{matrix of estimated (empirical Bayes) random effects, one row per group.}
#'     \item{lambda}{the penalty value used.}
#'     \item{loglik}{the observed-data log-likelihood at each iteration.}
#'     \item{AIC, BIC}{information criteria based on the number of non-zero coefficients.}
#'     \item{iter}{number of EM iterations performed.}
#'     \item{converged}{whether the EM algorithm reached `precision` and every
#'       M-step's internal `glmnet::glmnet()` lasso fit itself converged
#'       (see `glmnet_converged`). `FALSE` if either failed.}
#'     \item{glmnet_converged}{whether every M-step's `glmnet::glmnet()` call
#'       converged within its iteration cap. Coordinate descent can fail to
#'       converge for very small `lambda` combined with many covariates; when
#'       that happens here, the fit is stopped and the last iteration that
#'       did converge is returned instead of one corrupted by glmnet's
#'       fallback empty-model result, and a warning is issued. `FALSE` means
#'       this happened, and the fit at this `lambda` should not be trusted or
#'       compared to others without addressing it (e.g. a coarser lambda
#'       grid).}
#'     \item{x, z, y, nj, groups}{the design matrices/response/group sizes used to fit the model.}
#'     \item{call, formFixed, formRandom, groupVar}{the matched call and model specification.}
#'     \item{refit_fit}{if `refit = TRUE`, the `merMod` object from the
#'       post-selection `lme4::lmer()` refit (`NULL` if it failed -- see
#'       `refit`); absent entirely if `refit = FALSE`.}
#'     \item{refit_coefficients}{if `refit = TRUE`, named vector of refit
#'       fixed-effect point estimates, same names/order as `coefficients`,
#'       `NA` for covariates the lasso fit did not select (not `0` -- they
#'       were excluded from the refit model entirely, not estimated at
#'       zero). Access via `coef(fit, refit = TRUE)`.}
#'     \item{refit_singular}{if `refit = TRUE`, whether `lme4::isSingular()`
#'       flagged the refit as a boundary (singular) random-effects fit.}
#'   }
#'
#' @seealso [EMLasso_BIC_Search()] to select `lambda` automatically.
#'
#' @examples
#' set.seed(1)
#' n <- 30; nj1 <- 4; p <- 5; pstar <- 2
#' dat <- simulate_lasso_lmm(n = n, nj = nj1, p = p, pstar = pstar, seed = 1)
#' fixedForm <- stats::as.formula(paste("y ~", paste0("x", seq_len(p), collapse = " + ")))
#' fit <- EMLasso_LMM(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
#'                    lambda = 0.1, MaxIter = 50, verbose = FALSE)
#' fit
#' coef(fit)
#'
#' @export
EMLasso_LMM <- function(data, formFixed, groupVar, formRandom = ~1, lambda,
                         intercept = FALSE, precision = 1e-6, MaxIter = 200,
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

  # --- Standardize x (scale only, no centering) ---
  has_intercept <- "(Intercept)" %in% colnames(x)
  if (standardize_x) {
    x_sd <- apply(x, 2, stats::sd)
    if (has_intercept) x_sd[["(Intercept)"]] <- 1
    x_sd[x_sd == 0] <- 1 # guard against a constant/degenerate covariate
    x_orig <- x
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

    if (count == MaxIter) criterio <- precision / 10
  }

  teta <- c(beta, sigmae, D[upper.tri(D, diag = TRUE)])
  npar <- length(teta) - sum(beta == 0)
  AIC <- -2 * loglik + 2 * npar
  BIC <- -2 * loglik + log(N) * npar

  # --- Undo the internal standardization ---
  if (standardize_x) {
    beta <- beta / x_sd
    x <- x_orig
  }

  ranef <- t(b)
  rownames(ranef) <- dp$groups
  colnames(ranef) <- colnames(z)
  names(beta) <- colnames(x)

  obj <- list(
    coefficients = beta,
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
  class(obj) <- "EMLasso_LMM"

  if (isTRUE(refit)) {
    refit_res <- refit_lmm_impl(obj)
    obj$refit_fit <- refit_res$refit_fit
    obj$refit_coefficients <- refit_res$refit_coefficients
    obj$refit_formula <- refit_res$refit_formula
    obj$refit_singular <- refit_res$refit_singular
  }

  obj
}

#' @export
print.EMLasso_LMM <- function(x, digits = 4, ...) {
  cat("EM algorithm for lasso-penalized linear mixed models\n")
  cat("Call:\n")
  print(x$call)
  active <- x$coefficients[x$coefficients != 0]
  cat("\nlambda:", signif(x$lambda, digits),
      " | active covariates:", length(active), "of", length(x$coefficients), "\n")
  cat("\nFixed effects (non-zero only):\n")
  if (length(active) > 0) print(round(active, digits)) else cat("  (none)\n")
  cat("\nRandom effects covariance (D):\n")
  print(round(x$D, digits))
  cat("\nResidual variance:", round(x$sigma2, digits), "\n")
  cat("logLik:", round(utils::tail(x$loglik, 1), digits),
      " AIC:", round(x$AIC, digits), " BIC:", round(x$BIC, digits), "\n")
  cat("EM iterations:", x$iter, if (!x$converged) "(did not converge)" else "", "\n")
  if (!is.null(x$refit_coefficients)) {
    refit_active <- x$refit_coefficients[!is.na(x$refit_coefficients)]
    cat("\nPost-selection refit (unpenalized lme4::lmer(), point estimates only --\n")
    cat("do NOT use its SEs/p-values/CIs, they ignore the selection step):\n")
    print(round(refit_active, digits))
    if (isTRUE(x$refit_singular)) cat("(refit is singular; treat with extra caution)\n")
  }
  invisible(x)
}

#' @export
summary.EMLasso_LMM <- function(object, digits = 4, ...) {
  active <- object$coefficients[object$coefficients != 0]
  cat("EM algorithm for lasso-penalized linear mixed models\n")
  cat("Call:\n"); print(object$call)
  cat("\nNumber of groups:", object$n, " | total observations:", object$N, "\n")
  cat("lambda:", signif(object$lambda, digits), "\n\n")
  cat("Selected fixed effects (", length(active), "of", length(object$coefficients),
      "non-zero):\n")
  if (length(active) > 0) print(round(active, digits)) else cat("  (none)\n")
  cat("\nShrunk-to-zero covariates:",
      paste(names(object$coefficients)[object$coefficients == 0], collapse = ", "), "\n")
  cat("\nRandom effects covariance (D):\n"); print(round(object$D, digits))
  cat("Residual variance:", round(object$sigma2, digits), "\n")
  cat("\nlogLik:", round(utils::tail(object$loglik, 1), digits),
      " AIC:", round(object$AIC, digits), " BIC:", round(object$BIC, digits),
      " (", object$npar, "effective parameters )\n")
  if (!is.null(object$refit_coefficients)) {
    refit_active <- object$refit_coefficients[!is.na(object$refit_coefficients)]
    cat("\nPost-selection refit (unpenalized lme4::lmer() on the selected covariates):\n")
    print(round(refit_active, digits))
    cat("Point estimates only -- this refit's SEs/p-values/CIs do not account for\n")
    cat("the variable-selection step and must not be used for inference.\n")
    if (isTRUE(object$refit_singular)) cat("(refit is singular; treat with extra caution)\n")
  }
  invisible(object)
}

#' Extract coefficients from an `EMLasso_LMM` (or `EMLasso_BIC_Search`) fit
#'
#' @param object an `"EMLasso_LMM"` object.
#' @param refit if `TRUE`, return values from the post-selection
#'   `lme4::lmer()` refit instead of the lasso fit's (see the `refit`
#'   argument to [EMLasso_LMM()]; the fit/search must have been created with
#'   `refit = TRUE`) -- for `coef()`, its fixed-effect point estimates (`NA`
#'   for covariates the lasso did not select); for `fitted()`, its fitted
#'   values (fixed and random effects combined); for `residuals()`, `y`
#'   minus those fitted values. Default `FALSE`.
#' @param ... unused.
#' @export
coef.EMLasso_LMM <- function(object, refit = FALSE, ...) {
  if (isTRUE(refit)) {
    if (is.null(object$refit_coefficients)) {
      stop("This fit does not have a refit attached: call EMLasso_LMM() (or ",
           "EMLasso_BIC_Search()) with `refit = TRUE` first.", call. = FALSE)
    }
    return(object$refit_coefficients)
  }
  object$coefficients
}

#' @rdname coef.EMLasso_LMM
#' @export
fitted.EMLasso_LMM <- function(object, refit = FALSE, ...) {
  if (isTRUE(refit)) {
    if (is.null(object$refit_fit)) {
      stop("This fit does not have a refit attached: call EMLasso_LMM() (or ",
           "EMLasso_BIC_Search()) with `refit = TRUE` first.", call. = FALSE)
    }
    return(stats::fitted(object$refit_fit))
  }
  bmat <- t(object$ranef)
  nj <- object$nj; n <- length(nj)
  end_indices <- cumsum(nj)
  start_indices <- c(1, end_indices[-n] + 1)
  yhat <- numeric(sum(nj))
  for (j in seq_len(n)) {
    idx <- start_indices[j]:end_indices[j]
    yhat[idx] <- object$x[idx, , drop = FALSE] %*% object$coefficients +
      object$z[idx, , drop = FALSE] %*% bmat[, j]
  }
  yhat
}

#' @rdname coef.EMLasso_LMM
#' @export
residuals.EMLasso_LMM <- function(object, refit = FALSE, ...) {
  object$y - stats::fitted(object, refit = refit)
}
