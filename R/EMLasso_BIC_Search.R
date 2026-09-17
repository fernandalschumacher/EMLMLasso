#' Select the lasso penalty by BIC over a grid
#'
#' Internal engine, called by [EMLMLasso()] when it is asked to search over a
#' grid rather than fit a single lambda. Not exported -- use [EMLMLasso()].
#'
#' Fits the internal single-lambda engine over a grid of `lambda` values
#' (either a data-driven `glmnet`-generated grid, or one supplied directly
#' via `lambda_grid`) and selects the value that minimizes the Bayesian
#' Information Criterion (BIC).
#'
#' @param data,formFixed,groupVar,formRandom,intercept,precision,MaxIter,penalty_factor,standardize_x,refit,verbose
#'   passed straight through to the internal single-lambda engine at every
#'   grid point -- see [EMLMLasso()] for what each means.
#' @param n_lambdas number of `lambda` values in the data-driven grid (as in
#'   [glmnet::glmnet()]'s `nlambda`). Ignored if `lambda_grid` is supplied.
#' @param lambda_min_ratio smallest lambda in the data-driven grid, as a
#'   fraction of the largest (as in [glmnet::glmnet()]'s `lambda.min.ratio`).
#'   Ignored if `lambda_grid` is supplied.
#' @param lambda_grid optional numeric vector of lambda values to search
#'   directly, bypassing the data-driven grid construction (and
#'   `n_lambdas`/`lambda_min_ratio`) entirely. If `NULL` (the default), the
#'   grid is built the usual way, from an unpenalized-random-effects glmnet
#'   fit.
#'
#' @return a list with elements `best_model` (the winning internal single-
#'   lambda fit), `best_lambda`, `lambda_grid`, `bics`, `active_vars`,
#'   `converged` (per-lambda logical vector), `all_models`, `call`,
#'   `formFixed`, `formRandom`, `groupVar` -- see [EMLMLasso()] for what each
#'   means in the unified return object built from these.
#' @noRd
EMLasso_BIC_Search <- function(data, formFixed, groupVar, formRandom = ~1,
                                intercept = FALSE, n_lambdas = 30,
                                lambda_min_ratio = 0.001,
                                precision = 1e-4, MaxIter = 200,
                                penalty_factor = NULL, standardize_x = TRUE,
                                refit = FALSE, verbose = interactive(),
                                lambda_grid = NULL) {
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

  user_supplied_grid <- !is.null(lambda_grid)

  if (!user_supplied_grid) {
    # 1. Data-driven grid, from an unpenalized-random-effects glmnet fit
    has_intercept <- "(Intercept)" %in% colnames(x)
    init_x <- if (has_intercept) x[, colnames(x) != "(Intercept)", drop = FALSE] else x
    init_pf <- if (has_intercept) penalty_factor[colnames(x) != "(Intercept)"] else penalty_factor
    init_mod <- glmnet::glmnet(init_x, y, family = "gaussian", intercept = TRUE, alpha = 1,
                               nlambda = n_lambdas, lambda.min.ratio = lambda_min_ratio,
                               penalty.factor = init_pf)
    lambda_grid <- init_mod$lambda
  }

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
    if (user_supplied_grid) {
      warning("The BIC-optimal lambda is the smallest value in your `lambda` grid, ",
              "and BIC was still decreasing at that end -- the true optimum likely ",
              "lies at an even smaller lambda than this grid reached. Consider adding ",
              "smaller values to `lambda`.", call. = FALSE)
    } else {
      warning("The BIC-optimal lambda is the smallest value in the grid, and BIC was ",
              "still decreasing at that end -- the true optimum likely lies at an even ",
              "smaller lambda than this grid reached. Consider a smaller `lambda_min_ratio` ",
              "and/or a larger `n_lambdas`.", call. = FALSE)
    }
  }

  if (verbose) {
    cat("========================================\n")
    cat("Optimal model found at index", best_idx, "\n")
    cat("Best lambda:", lambda_grid[best_idx], "\n")
    cat("Active covariates:", active_vars[best_idx], "\n")
    cat("Lowest BIC:", bics[best_idx], "\n")
  }

  # 4. Refit only the selected model, not the whole grid (this is the
  # expensive, optional step -- see EMLMLasso()'s `refit` argument).
  if (isTRUE(refit)) {
    if (verbose) cat("\nRefitting the selected model (lme4::lmer(), unpenalized)...\n")
    refit_res <- refit_lmm_impl(best_model)
    best_model$refit_fit <- refit_res$refit_fit
    best_model$refit_coefficients <- refit_res$refit_coefficients
    best_model$refit_formula <- refit_res$refit_formula
    best_model$refit_singular <- refit_res$refit_singular
  }

  list(
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
}
