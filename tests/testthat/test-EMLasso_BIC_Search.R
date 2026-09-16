test_that("EMLasso_BIC_Search returns a consistent grid and selects its minimum", {
  dat <- simulate_lasso_lmm(n = 40, nj = 4, p = 6, pstar = 2, seed = 5)
  fixedForm <- as.formula(paste("y ~", paste0("x", 1:6, collapse = " + ")))
  search <- EMLasso_BIC_Search(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
                               n_lambdas = 6, MaxIter = 150, verbose = FALSE)

  expect_s3_class(search, "EMLasso_BIC_Search")
  expect_length(search$lambda_grid, 6)
  expect_length(search$bics, 6)
  expect_length(search$all_models, 6)
  expect_equal(search$best_lambda, search$lambda_grid[which.min(search$bics)])
  expect_s3_class(search$best_model, "EMLasso_LMM")
  expect_identical(coef(search), coef(search$best_model))
})

test_that("EMLasso_BIC_Search shrinks more at the largest grid lambda than at the smallest", {
  dat <- simulate_lasso_lmm(n = 30, nj = 4, p = 6, pstar = 2, seed = 8)
  fixedForm <- as.formula(paste("y ~", paste0("x", 1:6, collapse = " + ")))
  search <- EMLasso_BIC_Search(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
                               n_lambdas = 6, MaxIter = 150, verbose = FALSE)
  # weaker than strict path-monotonicity (not guaranteed here since each lambda is a
  # separate EM fit), but the smallest lambda should not shrink *more* than the largest.
  expect_true(search$active_vars[which.min(search$lambda_grid)] >=
                search$active_vars[which.max(search$lambda_grid)])
})

test_that("EMLasso_BIC_Search reports per-lambda convergence and excludes failures from BIC", {
  dat <- simulate_lasso_lmm(n = 30, nj = 4, p = 6, pstar = 2, seed = 8)
  fixedForm <- as.formula(paste("y ~", paste0("x", 1:6, collapse = " + ")))
  search <- EMLasso_BIC_Search(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
                               n_lambdas = 6, MaxIter = 150, verbose = FALSE)

  expect_length(search$converged, 6)
  expect_true(is.logical(search$converged))
  # a well-conditioned, low-dimensional example like this should converge everywhere
  expect_true(all(search$converged))
  expect_false(anyNA(search$bics))
})

test_that("a smaller lambda_min_ratio extends the grid to smaller lambda values", {
  dat <- simulate_lasso_lmm(n = 30, nj = 4, p = 6, pstar = 2, seed = 8)
  fixedForm <- as.formula(paste("y ~", paste0("x", 1:6, collapse = " + ")))
  search_default <- EMLasso_BIC_Search(dat, formFixed = fixedForm, groupVar = "id",
                                       formRandom = ~time, n_lambdas = 6, MaxIter = 100,
                                       verbose = FALSE)
  search_wider <- EMLasso_BIC_Search(dat, formFixed = fixedForm, groupVar = "id",
                                     formRandom = ~time, n_lambdas = 6, MaxIter = 100,
                                     lambda_min_ratio = 1e-5, verbose = FALSE)
  expect_true(min(search_wider$lambda_grid) < min(search_default$lambda_grid))
})

test_that("refit = TRUE refits only the selected model, not the whole grid", {
  skip_if_not_installed("lme4")
  dat <- simulate_lasso_lmm(n = 40, nj = 4, p = 6, pstar = 2, seed = 5)
  fixedForm <- as.formula(paste("y ~", paste0("x", 1:6, collapse = " + ")))
  search <- EMLasso_BIC_Search(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
                               n_lambdas = 6, MaxIter = 150, verbose = FALSE, refit = TRUE)

  expect_s4_class(search$best_model$refit_fit, "merMod")
  expect_equal(coef(search, refit = TRUE), search$best_model$refit_coefficients)
  # only the selected model was refit -- the rest of the grid was not touched
  other_models <- search$all_models[-which.min(search$bics)]
  expect_true(all(vapply(other_models, function(m) is.null(m$refit_fit), logical(1))))
})

test_that("formFixed = y ~ . works the same as in EMLasso_LMM()", {
  dat <- simulate_lasso_lmm(n = 30, nj = 4, p = 6, pstar = 2, seed = 8)
  search_dot <- EMLasso_BIC_Search(dat, formFixed = y ~ ., groupVar = "id", formRandom = ~time,
                                   n_lambdas = 6, MaxIter = 100, verbose = FALSE)
  search_explicit <- EMLasso_BIC_Search(
    dat, formFixed = y ~ time + x1 + x2 + x3 + x4 + x5 + x6, groupVar = "id",
    formRandom = ~time, n_lambdas = 6, MaxIter = 100, verbose = FALSE
  )
  expect_equal(deparse(search_dot$formFixed),
               "y ~ time + x1 + x2 + x3 + x4 + x5 + x6 - 1") # no intercept by default
  expect_equal(coef(search_dot), coef(search_explicit))
})

test_that("intercept is forwarded to every fit in the grid, default FALSE", {
  dat <- simulate_lasso_lmm(n = 30, nj = 4, p = 6, pstar = 2, seed = 8)
  fixedForm <- as.formula(paste("y ~", paste0("x", 1:6, collapse = " + ")))
  search_default <- EMLasso_BIC_Search(dat, formFixed = fixedForm, groupVar = "id",
                                       formRandom = ~time, n_lambdas = 4, MaxIter = 100,
                                       verbose = FALSE)
  search_intercept <- EMLasso_BIC_Search(dat, formFixed = fixedForm, groupVar = "id",
                                         formRandom = ~time, n_lambdas = 4, MaxIter = 100,
                                         intercept = TRUE, verbose = FALSE)
  expect_false("(Intercept)" %in% names(coef(search_default)))
  expect_true("(Intercept)" %in% names(coef(search_intercept)))
  expect_true(all(vapply(search_intercept$all_models,
                         function(m) "(Intercept)" %in% names(m$coefficients), logical(1))))
})
