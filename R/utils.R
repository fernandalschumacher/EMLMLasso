#' Build the response vector and design matrices used by the EM algorithm
#'
#' Internal helper that turns a `formFixed`/`formRandom`/`groupVar` + `data`
#' specification (the interface used throughout this package) into the
#' plain matrices the EM algorithm operates on. Rows are reordered (stably)
#' by `groupVar` so that observations from the same subject/cluster are
#' contiguous, which is what the algorithm assumes.
#'
#' @param formFixed two-sided formula for the fixed effects, e.g. `y ~ x1 + x2`.
#'   `y ~ .` is also accepted, meaning "every other column of `data` as a
#'   fixed-effect candidate", except the response and `groupVar` itself.
#' @param formRandom one-sided formula for the random effects, e.g. `~ 1` or `~ time`.
#' @param groupVar character string naming the subject/cluster variable in `data`.
#' @param data a data.frame containing all variables referenced above.
#' @param intercept if `TRUE`, force an intercept into the fixed-effects
#'   design matrix; if `FALSE` (the default), force it out. This always wins
#'   over whatever `formFixed` itself implies -- `y ~ x1 + x2` and
#'   `y ~ x1 + x2 - 1` are treated identically -- so it is the single place
#'   that controls fixed-effect intercept presence.
#'
#' @return a list with elements `y`, `x`, `z`, `nj` (observations per group,
#'   in the order groups first appear), `groups` (group labels in that same
#'   order), `order` (the row permutation applied to `data`), and the
#'   resolved `formFixed`/`formRandom` (with any `y ~ .` already expanded to
#'   the explicit covariate names it stands for, and the intercept forced in
#'   or out per `intercept`).
#' @noRd
prepare_lmm_data <- function(formFixed, formRandom, groupVar, data, intercept = FALSE) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  if (!is.character(groupVar) || length(groupVar) != 1L) {
    stop("`groupVar` must be a single character string naming a column of `data`.",
         call. = FALSE)
  }
  if (!groupVar %in% names(data)) {
    stop("`groupVar` = \"", groupVar, "\" was not found in `data`.", call. = FALSE)
  }
  formFixed <- stats::as.formula(formFixed)
  formRandom <- stats::as.formula(formRandom)
  if (length(formFixed) != 3L) {
    stop("`formFixed` must be a two-sided formula, e.g. `y ~ x1 + x2`.", call. = FALSE)
  }

  # `formFixed = y ~ .` requests every other column of `data` as a
  # fixed-effect candidate, except the response and `groupVar` itself.
  # (Any variable used only in `formRandom`, e.g. a random slope, is not
  # excluded -- it is swept in like any other column, matching how it would
  # be written out by hand if also wanted as a fixed effect.)
  if (identical(formFixed[[3]], quote(.))) {
    resp_vars <- all.vars(formFixed[[2]])
    candidates <- setdiff(names(data), c(resp_vars, groupVar))
    if (length(candidates) == 0) {
      stop("`formFixed` (`", deparse(formFixed), "`) expands to no fixed-effect ",
           "covariates: `data` has no columns left once the response and ",
           "`groupVar` are excluded.", call. = FALSE)
    }
    formFixed <- stats::as.formula(
      paste(deparse(formFixed[[2]]), "~", paste(candidates, collapse = " + "))
    )
  }

  # `intercept` always wins over whatever `formFixed` says: force it in (+1)
  # or out (+0/-1) via `update()` rather than post-processing the design
  # matrix, so factor contrasts (if any) are still handled correctly by
  # `model.matrix()` itself.
  formFixed <- stats::update(formFixed, if (isTRUE(intercept)) . ~ . + 1 else . ~ . + 0)

  vars_needed <- unique(c(all.vars(formFixed), all.vars(formRandom), groupVar))
  missing_vars <- setdiff(vars_needed, names(data))
  if (length(missing_vars) > 0) {
    stop("Variable(s) not found in `data`: ", paste(missing_vars, collapse = ", "),
         call. = FALSE)
  }
  data <- data[, vars_needed, drop = FALSE]
  if (anyNA(data)) {
    stop("`data` contains missing values in variables used by the model; ",
         "please remove or impute them before fitting.", call. = FALSE)
  }

  group_levels <- unique(data[[groupVar]])
  ord <- order(match(data[[groupVar]], group_levels))
  data <- data[ord, , drop = FALSE]

  mfFixed <- stats::model.frame(formFixed, data = data)
  y <- stats::model.response(mfFixed)
  x <- stats::model.matrix(formFixed, data = mfFixed)

  mfRandom <- stats::model.frame(formRandom, data = data)
  z <- stats::model.matrix(formRandom, data = mfRandom)

  nj <- as.integer(table(factor(data[[groupVar]], levels = group_levels)))

  list(
    y = as.numeric(y),
    x = x,
    z = z,
    nj = nj,
    groups = as.character(group_levels),
    order = ord,
    formFixed = formFixed,
    formRandom = formRandom
  )
}

#' Build a `thresh` argument compatible with both old and new glmnet
#'
#' glmnet >= 5.0 deprecated the top-level `thresh` (and `maxit`, `dfmax`,
#' `pmax`, `trace.it`) arguments in favor of `control = list(thresh = ...)`,
#' and emits a deprecation warning when the old form is used. Older glmnet
#' releases don't have a `control` formal at all. This picks whichever form
#' the installed glmnet supports, so neither version warns or misbehaves.
#' @noRd
glmnet_thresh_arg <- function(thresh) {
  if ("control" %in% names(formals(glmnet::glmnet))) {
    list(control = list(thresh = thresh))
  } else {
    list(thresh = thresh)
  }
}

#' Default penalty.factor: unpenalized intercept, penalized everything else
#' @noRd
default_penalty_factor <- function(x) {
  pf <- rep(1, ncol(x))
  if ("(Intercept)" %in% colnames(x)) {
    pf[colnames(x) == "(Intercept)"] <- 0
  }
  pf
}

#' Post-selection refit: an ordinary (unpenalized) `lme4::lmer()` fit on just
#' the lasso-selected covariates
#'
#' Lasso coefficients are shrunk toward zero by construction, which biases
#' them as point estimates even after the right covariates are selected.
#' This refits an ordinary linear mixed model -- no penalty -- on only the
#' covariates the lasso fit left non-zero (plus an intercept, if the
#' original `formFixed` had one), keeping the same random-effects structure.
#' It corrects the shrinkage bias in the fixed-effect point estimates, but
#' its standard errors/p-values/confidence intervals are **not** valid
#' inference: they treat the selected covariate set as fixed in advance,
#' when it was in fact chosen using the same data (see `refit` in
#' [EMLMLasso()] for the caveat this is meant to surface to callers).
#'
#' Building the refit data assumes `formFixed`/`formRandom` use simple,
#' untransformed variable names (as in every example in this package) --
#' factors, interactions or transformations in either formula are not
#' reconstructed from the stored design matrices and may not refit correctly.
#'
#' @param object a fitted model, as the list built by the internal `EMLasso_LMM()` engine.
#' @return a list with `refit_fit` (the fitted `merMod` object, or `NULL` if
#'   the refit failed), `refit_coefficients` (named numeric vector, same
#'   names/order as `object$coefficients`, `NA` for any covariate the lasso
#'   fit did not select), `refit_formula`, and `refit_singular` (whether
#'   `lme4::isSingular()` flagged a boundary fit).
#' @noRd
refit_lmm_impl <- function(object) {
  beta <- object$coefficients
  active <- names(beta)[beta != 0]
  has_intercept <- "(Intercept)" %in% colnames(object$x)
  fixed_terms <- setdiff(active, "(Intercept)")
  resp_name <- all.vars(object$formFixed)[1]

  xdf <- as.data.frame(object$x[, setdiff(colnames(object$x), "(Intercept)"), drop = FALSE])
  z_extra <- setdiff(setdiff(colnames(object$z), "(Intercept)"), colnames(xdf))
  if (length(z_extra) > 0) {
    xdf[z_extra] <- as.data.frame(object$z[, z_extra, drop = FALSE])
  }
  xdf[[object$groupVar]] <- factor(rep(object$groups, object$nj), levels = object$groups)
  xdf[[resp_name]] <- object$y

  rand_rhs <- deparse(object$formRandom[[2]])
  rand_term <- paste0("(", rand_rhs, " | ", object$groupVar, ")")
  fixed_rhs <- if (length(fixed_terms) == 0) "1" else paste(fixed_terms, collapse = " + ")
  if (!has_intercept && !isTRUE(object$centered)) fixed_rhs <- paste(fixed_rhs, "- 1")
  refit_formula <- stats::as.formula(paste(resp_name, "~", fixed_rhs, "+", rand_term))

  fit <- tryCatch(
    lme4::lmer(refit_formula, data = xdf, REML = TRUE,
               control = lme4::lmerControl(optimizer = "bobyqa", check.conv.singular = "ignore")),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    warning("Post-selection refit failed and was not attached to this fit: ",
            conditionMessage(fit), call. = FALSE)
    return(list(refit_fit = NULL, refit_coefficients = NULL,
                refit_formula = refit_formula, refit_singular = NA))
  }

  fe <- lme4::fixef(fit)
  refit_coefficients <- stats::setNames(rep(NA_real_, length(beta)), names(beta))
  common <- intersect(names(fe), names(refit_coefficients))
  refit_coefficients[common] <- fe[common]

  singular <- tryCatch(lme4::isSingular(fit), error = function(e) NA)
  if (isTRUE(singular)) {
    warning("Post-selection refit is singular (a boundary random-effects estimate); ",
            "point estimates may still be usable but treat with extra caution.",
            call. = FALSE)
  }

  list(refit_fit = fit, refit_coefficients = refit_coefficients,
       refit_formula = refit_formula, refit_singular = singular)
}

#' Observed-data log-likelihood for the linear mixed model
#'
#' `sum(log N(y_j; x_j beta, z_j D z_j' + sigmae I))` over clusters/subjects.
#' Used internally to monitor EM convergence and to compute AIC/BIC; it is
#' NOT penalized (the lasso penalty only enters the M-step).
#' @noRd
lmm_loglik <- function(beta, sigmae, D, y, x, z, nj) {
  n <- length(nj)
  p <- ncol(x)
  q <- ncol(z)
  end_indices <- cumsum(nj)
  start_indices <- c(1, end_indices[-n] + 1)
  suma1 <- 0
  for (j in seq_len(n)) {
    idx <- start_indices[j]:end_indices[j]
    y1 <- y[idx]
    x1 <- matrix(x[idx, ], ncol = p)
    z1 <- matrix(z[idx, ], ncol = q)
    med <- x1 %*% beta
    Psi <- z1 %*% D %*% t(z1) + sigmae * diag(nj[j])
    suma1 <- suma1 + mvtnorm::dmvnorm(y1, med, Psi, log = TRUE)
  }
  suma1
}
