test_that("EMLasso_LMM recovers a sparse fixed-effects structure", {
  dat <- simulate_lasso_lmm(n = 60, nj = 5, p = 8, pstar = 3, sigmae = 0.5, seed = 123)
  fixedForm <- as.formula(paste("y ~", paste0("x", 1:8, collapse = " + ")))
  fit <- EMLasso_LMM(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
                     lambda = 0.05, MaxIter = 300, precision = 1e-6, verbose = FALSE)

  expect_s3_class(fit, "EMLasso_LMM")
  expect_named(fit$coefficients, paste0("x", 1:8)) # no intercept by default
  # the true non-zero covariates should have the largest estimated coefficients
  true_active <- paste0("x", 1:3)
  other <- setdiff(paste0("x", 1:8), true_active)
  expect_true(min(abs(fit$coefficients[true_active])) >
                max(abs(fit$coefficients[other])))
  expect_true(fit$sigma2 > 0)
  expect_equal(dim(fit$D), c(2, 2))
  expect_true(fit$iter >= 1)
  expect_length(fitted(fit), nrow(dat))
  expect_length(residuals(fit), nrow(dat))
})

test_that("EMLasso_LMM leaves the intercept unpenalized when intercept = TRUE", {
  dat <- simulate_lasso_lmm(n = 20, nj = 4, p = 3, pstar = 1, seed = 7)
  fixedForm <- y ~ x1 + x2 + x3
  fit <- EMLasso_LMM(dat, formFixed = fixedForm, groupVar = "id", intercept = TRUE,
                     lambda = 5, MaxIter = 100, verbose = FALSE)
  intercept_idx <- which(names(fit$coefficients) == "(Intercept)")
  expect_equal(fit$penalty_factor[intercept_idx], 0)
})

test_that("EMLasso_LMM validates its inputs", {
  dat <- simulate_lasso_lmm(n = 10, nj = 3, p = 2, pstar = 1, seed = 3)
  expect_error(
    EMLasso_LMM(dat, formFixed = y ~ x1 + x2, groupVar = "not_a_column", lambda = 0.1),
    "not_a_column"
  )
  expect_error(
    EMLasso_LMM(dat, formFixed = y ~ x1 + x2, groupVar = "id", lambda = -1),
    "lambda"
  )
  expect_error(
    EMLasso_LMM(dat, formFixed = y ~ x1 + x2, groupVar = "id", lambda = 0.1,
               penalty_factor = c(1)),
    "penalty_factor"
  )
})

test_that("standardize_x leaves coefficients/fitted values essentially unchanged", {
  dat <- simulate_lasso_lmm(n = 40, nj = 4, p = 6, pstar = 2, sigmae = 0.5, seed = 21)
  fixedForm <- as.formula(paste("y ~", paste0("x", 1:6, collapse = " + ")))

  fit_std <- EMLasso_LMM(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
                         lambda = 0.1, MaxIter = 150, verbose = FALSE, standardize_x = TRUE)
  fit_raw <- EMLasso_LMM(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
                         lambda = 0.1, MaxIter = 150, verbose = FALSE, standardize_x = FALSE)

  # coefficients are reported back on the original scale of x either way
  expect_equal(fit_std$coefficients, fit_raw$coefficients, tolerance = 0.01)
  expect_equal(unname(fitted(fit_std)), unname(fitted(fit_raw)), tolerance = 0.01)
  # `x` stored on the fit is on the original (unstandardized) scale
  expect_equal(fit_std$x, fit_raw$x)
  expect_true(fit_std$glmnet_converged)
})

test_that("EMLasso_LMM does not require data pre-sorted by group", {
  set.seed(11)
  dat <- simulate_lasso_lmm(n = 15, nj = 3, p = 3, pstar = 1, seed = 11)
  shuffled <- dat[sample(nrow(dat)), ]

  fixedForm <- y ~ x1 + x2 + x3
  fit_sorted <- EMLasso_LMM(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
                            lambda = 0.2, MaxIter = 100, verbose = FALSE)
  fit_shuffled <- EMLasso_LMM(shuffled, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
                              lambda = 0.2, MaxIter = 100, verbose = FALSE)

  expect_equal(fit_sorted$coefficients, fit_shuffled$coefficients, tolerance = 1e-8)
})

test_that("refit attaches an unpenalized lme4 refit that reduces shrinkage bias", {
  skip_if_not_installed("lme4")
  dat <- simulate_lasso_lmm(n = 60, nj = 5, p = 8, pstar = 3, sigmae = 0.5, seed = 123)
  fixedForm <- as.formula(paste("y ~", paste0("x", 1:8, collapse = " + ")))
  fit <- EMLasso_LMM(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
                     lambda = 0.05, MaxIter = 300, precision = 1e-6, verbose = FALSE, refit = TRUE)

  expect_s4_class(fit$refit_fit, "merMod")
  expect_named(fit$refit_coefficients, names(fit$coefficients))
  active <- names(fit$coefficients)[fit$coefficients != 0]
  # the intercept is always in the refit (a location parameter, not treated
  # as "selected" or not), regardless of whether the lasso happened to zero
  # it out -- only the non-intercept, non-selected covariates should be NA
  inactive <- setdiff(names(fit$coefficients), union(active, "(Intercept)"))
  expect_false(anyNA(fit$refit_coefficients[active]))
  expect_true(all(is.na(fit$refit_coefficients[inactive])))
  # true non-zero coefficients are 1 in simulate_lasso_lmm(); the unpenalized
  # refit should land closer to that than the shrunk lasso estimate does
  true_active <- paste0("x", 1:3)
  expect_true(mean(abs(fit$refit_coefficients[true_active] - 1)) <
                mean(abs(fit$coefficients[true_active] - 1)))

  expect_equal(coef(fit, refit = TRUE), fit$refit_coefficients)
  expect_equal(coef(fit), fit$coefficients) # unaffected -- original kept alongside

  fit_norefit <- EMLasso_LMM(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
                             lambda = 0.05, MaxIter = 50, verbose = FALSE)
  expect_null(fit_norefit$refit_coefficients)
  expect_error(coef(fit_norefit, refit = TRUE), "refit")
})

test_that("refit handles a fully-sparse (0 active covariates) model", {
  skip_if_not_installed("lme4")
  dat <- simulate_lasso_lmm(n = 60, nj = 5, p = 8, pstar = 3, sigmae = 0.5, seed = 123)
  fixedForm <- as.formula(paste("y ~", paste0("x", 1:8, collapse = " + ")))
  fit0 <- EMLasso_LMM(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
                      lambda = 100, MaxIter = 50, verbose = FALSE, refit = TRUE,
                      intercept = TRUE)
  expect_equal(sum(fit0$coefficients != 0), 0)
  expect_s4_class(fit0$refit_fit, "merMod")
  expect_false(is.na(fit0$refit_coefficients[["(Intercept)"]]))
  expect_true(all(is.na(fit0$refit_coefficients[names(fit0$refit_coefficients) != "(Intercept)"])))
})

test_that("refit = TRUE errors without lme4 available", {
  dat <- simulate_lasso_lmm(n = 20, nj = 4, p = 3, pstar = 1, seed = 7)
  testthat::local_mocked_bindings(
    requireNamespace = function(...) FALSE, .package = "base"
  )
  expect_error(
    EMLasso_LMM(dat, formFixed = y ~ x1 + x2 + x3, groupVar = "id",
               lambda = 0.1, refit = TRUE, verbose = FALSE),
    "lme4"
  )
})

test_that("formFixed = y ~ . expands to every other column except the response and groupVar", {
  dat <- simulate_lasso_lmm(n = 20, nj = 4, p = 5, pstar = 2, seed = 7)
  fit_dot <- EMLasso_LMM(dat, formFixed = y ~ ., groupVar = "id", formRandom = ~time,
                         lambda = 0.1, MaxIter = 100, verbose = FALSE)
  fit_explicit <- EMLasso_LMM(dat, formFixed = y ~ time + x1 + x2 + x3 + x4 + x5,
                              groupVar = "id", formRandom = ~time,
                              lambda = 0.1, MaxIter = 100, verbose = FALSE)
  # `time` (used only in formRandom) is swept in like any other column,
  # `id` (groupVar) and `y` (the response) are excluded
  expect_equal(names(fit_dot$coefficients), names(fit_explicit$coefficients))
  expect_equal(fit_dot$coefficients, fit_explicit$coefficients)
  expect_equal(deparse(fit_dot$formFixed),
               "y ~ time + x1 + x2 + x3 + x4 + x5 - 1") # no intercept by default
})

test_that("formFixed = y ~ . errors when no covariates remain", {
  dat <- data.frame(y = rnorm(10), id = rep(1:2, 5))
  expect_error(
    EMLasso_LMM(dat, formFixed = y ~ ., groupVar = "id", lambda = 0.1, verbose = FALSE),
    "no fixed-effect"
  )
})

test_that("intercept defaults to FALSE, removing an implicit intercept", {
  dat <- simulate_lasso_lmm(n = 20, nj = 4, p = 3, pstar = 1, seed = 7)
  fixedForm <- y ~ x1 + x2 + x3 # implicit intercept in the formula itself
  fit <- EMLasso_LMM(dat, formFixed = fixedForm, groupVar = "id",
                     lambda = 0.1, MaxIter = 100, verbose = FALSE)
  expect_false("(Intercept)" %in% names(fit$coefficients))
  expect_false("(Intercept)" %in% colnames(fit$x))
})

test_that("intercept always wins over whatever formFixed says", {
  dat <- simulate_lasso_lmm(n = 20, nj = 4, p = 3, pstar = 1, seed = 7)
  fit_implicit <- EMLasso_LMM(dat, formFixed = y ~ x1 + x2 + x3, groupVar = "id",
                              lambda = 0.1, MaxIter = 100, intercept = TRUE, verbose = FALSE)
  fit_explicit_no <- EMLasso_LMM(dat, formFixed = y ~ x1 + x2 + x3 - 1, groupVar = "id",
                                 lambda = 0.1, MaxIter = 100, intercept = TRUE, verbose = FALSE)
  expect_true("(Intercept)" %in% names(fit_implicit$coefficients))
  expect_equal(fit_implicit$coefficients, fit_explicit_no$coefficients)

  fit_no_intercept_a <- EMLasso_LMM(dat, formFixed = y ~ x1 + x2 + x3, groupVar = "id",
                                    lambda = 0.1, MaxIter = 100, intercept = FALSE, verbose = FALSE)
  fit_no_intercept_b <- EMLasso_LMM(dat, formFixed = y ~ x1 + x2 + x3 - 1, groupVar = "id",
                                    lambda = 0.1, MaxIter = 100, intercept = FALSE, verbose = FALSE)
  expect_false("(Intercept)" %in% names(fit_no_intercept_a$coefficients))
  expect_equal(fit_no_intercept_a$coefficients, fit_no_intercept_b$coefficients)
})

test_that("fitted/residuals(refit = TRUE) use the lme4 refit's own values", {
  skip_if_not_installed("lme4")
  dat <- simulate_lasso_lmm(n = 60, nj = 5, p = 8, pstar = 3, sigmae = 0.5, seed = 123)
  fixedForm <- as.formula(paste("y ~", paste0("x", 1:8, collapse = " + ")))
  fit <- EMLasso_LMM(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
                     lambda = 0.05, MaxIter = 300, precision = 1e-6, verbose = FALSE, refit = TRUE)

  expect_equal(fitted(fit, refit = TRUE), stats::fitted(fit$refit_fit), ignore_attr = TRUE)
  expect_equal(residuals(fit, refit = TRUE), fit$y - fitted(fit, refit = TRUE), ignore_attr = TRUE)
  # the two should generally differ from the lasso fit's own fitted values
  expect_false(isTRUE(all.equal(fitted(fit, refit = TRUE), fitted(fit))))

  fit_norefit <- EMLasso_LMM(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
                             lambda = 0.05, MaxIter = 50, verbose = FALSE)
  expect_error(fitted(fit_norefit, refit = TRUE), "refit")
  expect_error(residuals(fit_norefit, refit = TRUE), "refit")
})
