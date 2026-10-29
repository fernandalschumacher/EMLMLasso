#' Fit a lasso-penalized linear mixed model via EM
#'
#' Fits a linear mixed model with a lasso penalty on the fixed-effects
#' coefficients, using an Expectation-Maximization (EM) algorithm (at each
#' E-step, random effects are integrated out in closed form; at each M-step,
#' the fixed effects are updated by a coordinate-descent lasso fit
#' ([glmnet::glmnet()]) on the random-effects-adjusted response, and the
#' variance components are updated in closed form). This is the package's main 
#' function and the argument `lambda` decides whether it fits one model
#' or searches a grid and selects by BIC.
#'
#' Following Oliveira et al. (2026), the model is
#' \eqn{Y_i = X_i \beta + Z_i b_i + \epsilon_i}, with
#' \eqn{b_i \sim N_q(0, D)} independent of \eqn{\epsilon_i \sim N_{n_i}(0,
#' \sigma^2 I_{n_i})}, so the covariance for subject/cluster
#' \eqn{i} is \eqn{\Psi_i = Z_i D Z_i' + \sigma^2 I_{n_i}}. `lambda` plays the
#' role of \eqn{\lambda^{(k)}} in the M-step of Algorithm 1 of Oliveira et al.
#' (2026) (their Equation 8) -- the penalty passed directly to
#' [glmnet::glmnet()]. Unlike Algorithm 1, a `lambda` value is fixed
#' for the duration of one EM fit rather than re-selected by
#' cross-validation at every EM iteration; see `NEWS.md` for why this
#' package selects it by a BIC grid search instead. Following Equation 1 of
#' the same paper, no intercept is fit by default (see `intercept`).
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
#' @param lambda controls what gets fit:
#'   \describe{
#'     \item{`NULL` (the default)}{search a data-driven grid of `n_lambdas`
#'       values (as in [glmnet::glmnet()]'s own grid construction, spanning
#'       down to `lambda_min_ratio` of the largest value) and select by BIC.}
#'     \item{a single non-negative number}{fit that one lambda directly, no
#'       search, no BIC comparison.}
#'     \item{a numeric vector of two or more non-negative values}{search
#'       exactly this grid (sorted, deduplicated, descending) and select by
#'       BIC, instead of building one automatically -- `n_lambdas` and
#'       `lambda_min_ratio` are ignored in this case.}
#'   }
#' @param intercept if `TRUE`, force an intercept into the fixed-effects
#'   design matrix; if `FALSE` (the default), force it out. This always wins
#'   over whatever `formFixed` itself implies -- `y ~ x1 + x2` and
#'   `y ~ x1 + x2 - 1` are treated identically.
#' @param n_lambdas number of `lambda` values in the data-driven grid (as in
#'   [glmnet::glmnet()]'s `nlambda`). Only used when `lambda = NULL`. Default `30`.
#' @param lambda_min_ratio smallest lambda in the data-driven grid, as a
#'   fraction of the largest (as in [glmnet::glmnet()]'s `lambda.min.ratio`).
#'   Only used when `lambda = NULL`. Default `0.001`, ten times smaller than
#'   glmnet's own default for `p > n` data (`0.01`).
#' @param precision convergence tolerance on the Euclidean distance between
#'   successive parameter vectors (beta, sigma2, D). Default `1e-4`.
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
#'   `fitted()`, etc. are unaffected by this option. Set to `FALSE` to fit on
#'   the original scale of `x` instead (relying only on glmnet's own default
#'   internal standardization, recomputed every M-step).
#' @param refit if `TRUE`, once the (selected or only) model is available,
#'   also refits an ordinary (unpenalized) `lme4::lmer()` model on just its
#'   non-zero covariates -- same random-effects structure -- to correct the
#'   shrinkage bias lasso point estimates always carry. When `lambda`
#'   triggers a search, only the selected (best) model is refit, not every
#'   grid point. Requires the `lme4` package. **The refit's own standard
#'   errors, p-values and confidence intervals must not be used or
#'   reported** -- the covariate set was chosen using this same data, so
#'   they do not account for that selection step and understate the true
#'   uncertainty (classic post-selection inference). The original lasso fit
#'   (`coefficients`, `AIC`, `BIC`, etc.) is always kept unchanged alongside
#'   it. Default `FALSE`.
#' @param verbose if `TRUE`, prints progress (one line every 10 EM
#'   iterations, plus, when searching, one line per grid value). Defaults to
#'   `interactive()`.
#'
#' @return An object of class `"EMLMLasso"`, a list with elements:
#'   \describe{
#'     \item{coefficients}{named vector of estimated fixed-effects
#'       coefficients for the selected (or only) model.}
#'     \item{sigma2}{estimated residual variance.}
#'     \item{D}{estimated random-effects covariance matrix.}
#'     \item{ranef}{matrix of estimated (empirical Bayes) random effects, one row per group.}
#'     \item{lambda}{the penalty value of the selected (or only) model.}
#'     \item{loglik}{the observed-data log-likelihood at each iteration, for
#'       the selected (or only) model.}
#'     \item{AIC, BIC}{information criteria for the selected (or only) model.}
#'     \item{iter}{number of EM iterations performed for the selected (or only) model.}
#'     \item{converged}{whether the selected (or only) model's EM algorithm
#'       reached `precision` and its M-steps' `glmnet::glmnet()` fits
#'       themselves converged (see `glmnet_converged`).}
#'     \item{glmnet_converged}{whether every M-step of the selected (or
#'       only) model's `glmnet::glmnet()` calls converged.}
#'     \item{x, z, y, nj, groups}{the design matrices/response/group sizes used to fit the model.}
#'     \item{call, formFixed, formRandom, groupVar}{the matched call and model specification.}
#'     \item{search_performed}{`TRUE` if `lambda` triggered a grid search
#'       (`NULL` or a vector), `FALSE` if a single lambda was fit directly.}
#'     \item{lambda_grid}{the full grid of lambda values considered, or
#'       `NULL` if `search_performed` is `FALSE`.}
#'     \item{bics, active_vars}{BIC and number of active (non-zero)
#'       covariates at each grid value (`NA` at any lambda excluded because
#'       its M-step's `glmnet::glmnet()` fit did not converge), or `NULL` if
#'       `search_performed` is `FALSE`.}
#'     \item{grid_converged}{logical vector, one per grid value: whether
#'       that lambda's fit converged and was used in the BIC comparison, or
#'       `NULL` if `search_performed` is `FALSE`.}
#'     \item{all_models}{list of all fits considered, one per grid value
#'       (including any excluded from BIC comparison, for inspection), or
#'       `NULL` if `search_performed` is `FALSE`.}
#'   }
#'   If `refit = TRUE`, this additionally carries `refit_fit`,
#'   `refit_coefficients`, `refit_formula` and `refit_singular` for the
#'   selected (or only) model -- see the `refit` argument above.
#'
#' @details
#' The BIC grid search (`lambda = NULL` or a vector) replaces an earlier
#' cross-validation-based selection (`cv.glmnet` on the random-effects-
#' adjusted response, refit at every EM iteration -- the `While` loop's
#' "Obtain the optimum value \eqn{\lambda^{(k)}} via cross-validation" step
#' in Algorithm 1 of Oliveira et al. 2026), which was both much slower and,
#' in simulations, noticeably less stable than selecting by BIC on a fixed
#' grid. See `NEWS.md` for details.
#'
#' @examples
#' set.seed(1)
#' p <- 5
#' dat <- simulate_lasso_lmm(n = 25, nj = 4, p = p, pstar = 2, seed = 1)
#' fixedForm <- stats::as.formula(paste("y ~", paste0("x", seq_len(p), collapse = " + ")))
#'
#' # fit a single lambda directly
#' fit <- EMLMLasso(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
#'                  lambda = 0.1, MaxIter = 50, verbose = FALSE)
#' fit
#' coef(fit)
#'
#' # search a data-driven grid and select by BIC (the common case)
#' search <- EMLMLasso(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
#'                     n_lambdas = 8, MaxIter = 50, verbose = FALSE)
#' search
#' coef(search)
#'
#' @export
EMLMLasso <- function(data, formFixed, groupVar, formRandom = ~1, lambda = NULL,
                       intercept = FALSE, n_lambdas = 30, lambda_min_ratio = 0.001,
                       precision = 1e-4, MaxIter = 200, penalty_factor = NULL,
                       standardize_x = TRUE, refit = FALSE, verbose = interactive()) {
  cl <- match.call()
  if (!is.null(lambda)) {
    if (!is.numeric(lambda) || length(lambda) == 0L || anyNA(lambda) || any(lambda < 0)) {
      stop("`lambda` must be NULL, a single non-negative number, or a numeric ",
           "vector of non-negative values.", call. = FALSE)
    }
  }
  if (isTRUE(refit) && !requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for `refit = TRUE`. Install it with ",
         "install.packages(\"lme4\").", call. = FALSE)
  }

  search_performed <- is.null(lambda) || length(lambda) > 1L

  if (!search_performed) {
    fit <- EMLasso_LMM(data = data, formFixed = formFixed, groupVar = groupVar,
                        formRandom = formRandom, lambda = lambda, intercept = intercept,
                        precision = precision, MaxIter = MaxIter,
                        penalty_factor = penalty_factor, standardize_x = standardize_x,
                        refit = refit, verbose = verbose)
    obj <- c(fit, list(
      search_performed = FALSE, lambda_grid = NULL, bics = NULL,
      active_vars = NULL, grid_converged = NULL, all_models = NULL
    ))
  } else {
    user_grid <- if (!is.null(lambda)) sort(unique(lambda), decreasing = TRUE) else NULL
    search <- EMLasso_BIC_Search(data = data, formFixed = formFixed, groupVar = groupVar,
                                  formRandom = formRandom, intercept = intercept,
                                  n_lambdas = n_lambdas, lambda_min_ratio = lambda_min_ratio,
                                  precision = precision, MaxIter = MaxIter,
                                  penalty_factor = penalty_factor, standardize_x = standardize_x,
                                  refit = refit, verbose = verbose, lambda_grid = user_grid)
    obj <- c(search$best_model, list(
      search_performed = TRUE, lambda_grid = search$lambda_grid, bics = search$bics,
      active_vars = search$active_vars, grid_converged = search$converged,
      all_models = search$all_models
    ))
  }
  obj$call <- cl
  class(obj) <- "EMLMLasso"
  obj
}

#' @export
print.EMLMLasso <- function(x, digits = 4, ...) {
  cat("EM algorithm for lasso-penalized linear mixed models\n")
  cat("Call:\n"); print(x$call)
  if (isTRUE(x$search_performed)) {
    cat("\nBIC search over", length(x$lambda_grid), "lambda values\n")
    if (!is.null(x$grid_converged) && any(!x$grid_converged)) {
      cat(sum(!x$grid_converged), "of", length(x$lambda_grid), "lambda value(s) excluded ",
          "from comparison (glmnet did not converge; see `grid_converged`).\n")
    }
  }
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
summary.EMLMLasso <- function(object, digits = 4, ...) {
  active <- object$coefficients[object$coefficients != 0]
  cat("EM algorithm for lasso-penalized linear mixed models\n")
  cat("Call:\n"); print(object$call)
  cat("\nNumber of groups:", object$n, " | total observations:", object$N, "\n")
  if (isTRUE(object$search_performed)) {
    cat("BIC search over", length(object$lambda_grid), "lambda values\n")
    #cat("lambda grid: ", signif(min(object$lambda_grid), digits), " to ",
    #    signif(max(object$lambda_grid), digits), " (", length(object$lambda_grid),
    #    " values)\n", sep = "")
    #cat("Active covariates by lambda:", paste(object$active_vars, collapse = ", "), "\n\n")
    cat("\nSelected lambda (by BIC):", signif(object$lambda, digits), "\n\n")
  } else {
    cat("lambda:", signif(object$lambda, digits), "\n\n")
  }
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

#' Extract coefficients, fitted values or residuals from an `EMLMLasso` fit
#'
#' @param object an `"EMLMLasso"` object.
#' @param refit if `TRUE`, return values from the post-selection
#'   `lme4::lmer()` refit instead of the lasso fit's (see the `refit`
#'   argument to [EMLMLasso()]; the object must have been created with
#'   `refit = TRUE`) -- for `coef()`, its fixed-effect point estimates (`NA`
#'   for covariates the lasso did not select); for `fitted()`, its fitted
#'   values (fixed and random effects combined); for `residuals()`, `y`
#'   minus those fitted values. Default `FALSE`.
#' @param ... unused.
#' @return For `coef()`, a named numeric vector of fixed-effect estimates,
#'   one per column of the fixed-effects design matrix (including the
#'   intercept when `intercept = TRUE` in [EMLMLasso()]); covariates shrunk
#'   out by the lasso are exactly `0`, or `NA` when `refit = TRUE`.
#'   For `fitted()` and `residuals()`, a numeric vector of length equal to
#'   the total number of observations, in the same order as the stored
#'   response `y` (observations grouped by `groupVar`), with fitted values
#'   combining the fixed effects and the predicted random effects.
#' @export
coef.EMLMLasso <- function(object, refit = FALSE, ...) {
  if (isTRUE(refit)) {
    if (is.null(object$refit_coefficients)) {
      stop("This fit does not have a refit attached: call EMLMLasso() with ",
           "`refit = TRUE` first.", call. = FALSE)
    }
    return(object$refit_coefficients)
  }
  object$coefficients
}

#' @rdname coef.EMLMLasso
#' @export
fitted.EMLMLasso <- function(object, refit = FALSE, ...) {
  if (isTRUE(refit)) {
    if (is.null(object$refit_fit)) {
      stop("This fit does not have a refit attached: call EMLMLasso() with ",
           "`refit = TRUE` first.", call. = FALSE)
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

#' @rdname coef.EMLMLasso
#' @export
residuals.EMLMLasso <- function(object, refit = FALSE, ...) {
  object$y - stats::fitted(object, refit = refit)
}

#' Plot BIC against the lambda grid
#'
#' Only applicable to an `"EMLMLasso"` object created with a BIC search
#' (`lambda = NULL` or a vector of candidate values, the default) -- errors
#' informatively if `x` was fit at a single `lambda` instead, since there is
#' then no grid to plot.
#'
#' @param x an `"EMLMLasso"` object.
#' @param ... further arguments passed to [graphics::plot()].
#' @return Invisibly returns `x`, the input `"EMLMLasso"` object, unchanged.
#'   Called for its side effect of plotting BIC against the lambda grid, with
#'   the selected lambda highlighted in red.
#' @export
plot.EMLMLasso <- function(x, ...) {
  if (!isTRUE(x$search_performed)) {
    stop("Nothing to plot: this object was fit at a single `lambda`, not a ",
         "BIC search over a grid. Call EMLMLasso() with `lambda = NULL` or a ",
         "vector of candidate values to get a plottable BIC search.", call. = FALSE)
  }
  best_idx <- which.min(x$bics)
  graphics::plot(x$lambda_grid, x$bics, type = "b", xlab = "lambda", ylab = "BIC", ...)
  graphics::points(x$lambda_grid[best_idx], x$bics[best_idx], pch = 19, col = "red")
  invisible(x)
}
