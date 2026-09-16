#' Select the lasso penalty by BIC over a data-driven grid
#'
#' Fits [EMLasso_LMM()] over a grid of `lambda` values (the same
#' `glmnet`-generated grid used to initialize each fit) and selects the
#' value that minimizes the Bayesian Information Criterion (BIC).
#'
#' @references
#' Oliveira DCR, Schumacher FL, Lachos VH (2026). The use of the EM algorithm
#' for regularization problems in high-dimensional linear mixed-effects
#' models. *Statistical Methods in Medical Research*, 35(2), 240-254.
#' \doi{10.1177/09622802251399913}
#'
#' @inheritParams EMLasso_LMM
#' @param n_lambdas number of `lambda` values in the grid (as in
#'   [glmnet::glmnet()]'s `nlambda`). Default `30`.
#' @param lambda_min_ratio smallest lambda in the grid, as a fraction of the
#'   largest (as in [glmnet::glmnet()]'s `lambda.min.ratio`). Default
#'   `0.001`, ten times smaller than glmnet's own default for `p > n` data
#'   (`0.01`).
#' @param refit if `TRUE`, once the best lambda is selected, also refits an
#'   ordinary (unpenalized) `lme4::lmer()` model on just its selected
#'   covariates -- see [EMLasso_LMM()]'s `refit` argument. Only the selected
#'   (best) model is refit, not the whole grid. Default `FALSE`.
#'
#' @return An object of class `"EMLasso_BIC_Search"`, a list with elements:
#'   \describe{
#'     \item{best_model}{the `"EMLasso_LMM"` fit at the selected lambda.}
#'     \item{best_lambda}{the selected penalty value.}
#'     \item{lambda_grid}{the full grid of lambda values considered.}
#'     \item{bics, active_vars}{BIC and number of active (non-zero) covariates
#'       at each grid value; `NA` at any lambda excluded from comparison
#'       because its M-step's `glmnet::glmnet()` fit did not converge (see
#'       `converged` and [EMLasso_LMM()]'s `glmnet_converged`).}
#'     \item{converged}{logical vector, one per grid value: whether that
#'       lambda's fit converged and was used in the BIC comparison.}
#'     \item{all_models}{list of all `"EMLasso_LMM"` fits, one per grid value
#'       (including any excluded from BIC comparison, for inspection).}
#'     \item{call, formFixed, formRandom, groupVar}{the matched call and model specification.}
#'   }
#'   If `refit = TRUE`, `best_model` additionally carries `refit_fit`,
#'   `refit_coefficients` and `refit_singular` -- see [EMLasso_LMM()]'s
#'   `refit` argument. Only the selected (best) model is refit, not the whole grid.
#'
#' @seealso [EMLasso_LMM()] to fit a single lambda directly.
#'
#' @details
#' This version replaces an earlier cross-validation-based selection
#' (`cv.glmnet` on the random-effects-adjusted response, refit at every
#' EM iteration -- the `While` loop's "Obtain the optimum value \eqn{\lambda^{(k)}} via
#' cross-validation" step in Algorithm 1 of Oliveira et al. 2026), which was
#' both much slower and, in simulations, noticeably less stable than
#' selecting by BIC on a fixed grid. See `NEWS.md` for details.
#'
#'
#' @examples
#' set.seed(1)
#' p <- 5
#' dat <- simulate_lasso_lmm(n = 25, nj = 4, p = p, pstar = 2, seed = 1)
#' fixedForm <- stats::as.formula(paste("y ~", paste0("x", seq_len(p), collapse = " + ")))
#' search <- EMLasso_BIC_Search(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
#'                              n_lambdas = 8, MaxIter = 50, verbose = FALSE)
#' search
#' coef(search)
#'
#' @export
EMLasso_BIC_Search <- function(data, formFixed, groupVar, formRandom = ~1,
                                intercept = FALSE, n_lambdas = 30,
                                lambda_min_ratio = 0.001,
                                precision = 1e-4, MaxIter = 200,
                                penalty_factor = NULL, standardize_x = TRUE,
                                refit = FALSE, verbose = interactive()) {
  cl <- match.call()
  if (isTRUE(refit) && !requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for `refit = TRUE`. Install it with ",
         "install.packages(\"lme4\").", call. = FALSE)
  }
  dp <- prepare_lmm_data(formFixed, formRandom, groupVar, data, intercept = intercept)
  y <- dp$y; x <- dp$x

  if (is.null(penalty_factor)) {
    penalty_factor <- default_penalty_factor(x)
  } else if (length(penalty_factor) != ncol(x)) {
    stop("`penalty_factor` must have length ", ncol(x), " (one entry per fixed-effect ",
         "coefficient, including the intercept if present).", call. = FALSE)
  }

  # 1. Data-driven grid, from an unpenalized-random-effects glmnet fit
  has_intercept <- "(Intercept)" %in% colnames(x)
  init_x <- if (has_intercept) x[, colnames(x) != "(Intercept)", drop = FALSE] else x
  init_pf <- if (has_intercept) penalty_factor[colnames(x) != "(Intercept)"] else penalty_factor
  init_mod <- glmnet::glmnet(init_x, y, family = "gaussian", intercept = TRUE, alpha = 1,
                             nlambda = n_lambdas, lambda.min.ratio = lambda_min_ratio,
                             penalty.factor = init_pf)
  lambda_grid <- init_mod$lambda

  if (verbose) {
    cat("Generated", length(lambda_grid), "data-driven lambdas ranging from",
        round(max(lambda_grid), 4), "to", round(min(lambda_grid), 4), "\n\n")
  }

  results <- vector("list", length(lambda_grid))
  bics <- numeric(length(lambda_grid))
  active_vars <- numeric(length(lambda_grid))
  converged <- logical(length(lambda_grid))

  # 2. Fit the EM algorithm over the grid
  for (i in seq_along(lambda_grid)) {
    current_lambda <- lambda_grid[i]
    if (verbose) cat("Fitting model", i, "| lambda:", round(current_lambda, 5), "\n")

    fit <- EMLasso_LMM(data = data, formFixed = formFixed, groupVar = groupVar,
                       formRandom = formRandom, lambda = current_lambda,
                       intercept = intercept, precision = precision, MaxIter = MaxIter,
                       penalty_factor = penalty_factor, standardize_x = standardize_x,
                       verbose = verbose)

    results[[i]] <- fit
    converged[i] <- fit$glmnet_converged
    bics[i] <- if (converged[i]) fit$BIC else NA_real_
    active_vars[i] <- if (converged[i]) sum(fit$coefficients != 0) else NA_real_

    if (verbose) {
      if (converged[i]) {
        cat("--> active covariates:", active_vars[i], "| BIC:", round(bics[i], 2), "\n\n")
      } else {
        cat("--> glmnet did not converge at this lambda; excluded from BIC comparison\n\n")
      }
    }
  }

  if (all(is.na(bics))) {
    stop("glmnet did not converge at any lambda in the grid; none of these fits ",
         "can be trusted. Try a coarser (larger) lambda grid, e.g. a smaller ",
         "`n_lambdas` or a larger `lambda.min.ratio` via a custom grid.",
         call. = FALSE)
  }
  if (any(!converged) && verbose) {
    cat("Note:", sum(!converged), "of", length(lambda_grid),
        "lambda value(s) were excluded from BIC comparison (glmnet did not converge).\n\n")
  }

  # 3. Select the best model among those that converged
  best_idx <- which.min(bics)
  best_model <- results[[best_idx]]

  if (best_idx == length(lambda_grid) && length(lambda_grid) > 1 && !is.na(bics[best_idx - 1]) &&
      bics[best_idx] < bics[best_idx - 1]) {
    warning("The BIC-optimal lambda is the smallest value in the grid, and BIC was ",
            "still decreasing at that end -- the true optimum likely lies at an even ",
            "smaller lambda than this grid reached. Consider a smaller `lambda_min_ratio` ",
            "and/or a larger `n_lambdas`.", call. = FALSE)
  }

  if (verbose) {
    cat("========================================\n")
    cat("Optimal model found at index", best_idx, "\n")
    cat("Best lambda:", lambda_grid[best_idx], "\n")
    cat("Active covariates:", active_vars[best_idx], "\n")
    cat("Lowest BIC:", bics[best_idx], "\n")
  }

  # 4. Refit only the selected model, not the whole grid (this is the
  # expensive, optional step -- see EMLasso_LMM()'s `refit` argument).
  if (isTRUE(refit)) {
    if (verbose) cat("\nRefitting the selected model (lme4::lmer(), unpenalized)...\n")
    refit_res <- refit_lmm_impl(best_model)
    best_model$refit_fit <- refit_res$refit_fit
    best_model$refit_coefficients <- refit_res$refit_coefficients
    best_model$refit_formula <- refit_res$refit_formula
    best_model$refit_singular <- refit_res$refit_singular
  }

  obj <- list(
    best_model = best_model,
    best_lambda = lambda_grid[best_idx],
    lambda_grid = lambda_grid,
    bics = bics,
    active_vars = active_vars,
    converged = converged,
    all_models = results,
    call = cl, formFixed = dp$formFixed,
    formRandom = dp$formRandom, groupVar = groupVar
  )
  class(obj) <- "EMLasso_BIC_Search"
  obj
}

#' @export
print.EMLasso_BIC_Search <- function(x, digits = 4, ...) {
  best_idx <- which.min(x$bics)
  cat("EM lasso BIC search over", length(x$lambda_grid), "lambda values\n")
  cat("Call:\n"); print(x$call)
  if (!is.null(x$converged) && any(!x$converged)) {
    cat("\n", sum(!x$converged), " of ", length(x$lambda_grid), " lambda value(s) excluded ",
        "from comparison (glmnet did not converge; see `converged`).\n", sep = "")
  }
  cat("\nBest lambda:", signif(x$best_lambda, digits),
      "| active covariates:", x$active_vars[best_idx],
      "| BIC:", round(x$bics[best_idx], digits), "\n\n")
  cat("Selected model:\n")
  print(x$best_model)
  invisible(x)
}

#' @export
summary.EMLasso_BIC_Search <- function(object, digits = 4, ...) {
  cat("EM lasso BIC search over", length(object$lambda_grid), "lambda values\n")
  cat("Call:\n"); print(object$call)
  #cat("\nlambda grid: ", signif(min(object$lambda_grid), digits), " to ",
  #    signif(max(object$lambda_grid), digits), " (", length(object$lambda_grid),
  #    " values)\n", sep = "")
  #cat("Active covariates by lambda:", paste(object$active_vars, collapse = ", "), "\n\n")
  cat("\nSelected model (by BIC):\n")
  summary(object$best_model, digits = digits)
  invisible(object)
}

#' @rdname coef.EMLasso_LMM
#' @export
coef.EMLasso_BIC_Search <- function(object, refit = FALSE, ...) {
  stats::coef(object$best_model, refit = refit)
}

#' @rdname coef.EMLasso_LMM
#' @export
fitted.EMLasso_BIC_Search <- function(object, refit = FALSE, ...) {
  stats::fitted(object$best_model, refit = refit)
}

#' Plot BIC against the lambda grid
#'
#' @param x an `"EMLasso_BIC_Search"` object.
#' @param ... further arguments passed to [graphics::plot()].
#' @export
plot.EMLasso_BIC_Search <- function(x, ...) {
  best_idx <- which.min(x$bics)
  graphics::plot(x$lambda_grid, x$bics, type = "b", xlab = "lambda", ylab = "BIC", ...)
  graphics::points(x$lambda_grid[best_idx], x$bics[best_idx], pch = 19, col = "red")
  invisible(x)
}
