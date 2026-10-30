# EMLMLasso 0.1.0

## CRAN resubmission

* New `center` argument to `EMLMLasso()` (default `TRUE`): with
  `intercept = FALSE` (the default), the response and covariates are now
  centred before fitting, matching the centred setting assumed by Oliveira
  et al. (2026). Coefficients stay on the original scale; the implied
  intercept is returned as `implied_intercept` and used by `fitted()` and
  `residuals()`, and the `lme4` refit includes an intercept. Use
  `center = FALSE` for the previous behaviour.
* Bug fix: a fit that stopped because it reached `MaxIter` was reported as
  `converged = TRUE`; it is now `FALSE`, and `print()` shows
  "(did not converge)".
* Documentation now states that this package modifies the published
  algorithm (lambda selected by BIC over a grid instead of per-iteration
  cross-validation; different convergence criterion), and aligns notation
  with the paper.
* The riboflavin vignette now standardizes the response and all
  covariates, as in Section 5 of the paper, and uses `MaxIter = 500` so
  that the selected fit converges.

* Documented the return values of `coef()`/`fitted()`/`residuals()` and
  `plot()` methods for `"EMLMLasso"` objects.
* Removed single quotes around function names in DESCRIPTION.

## Combining functions as `EMLMLasso()`

* `EMLasso_LMM()` and `EMLasso_BIC_Search()` are replaced by a single
   wrapper function, `EMLMLasso()`, matching the package name. What
  `lambda` is now decides the mode: `NULL` (the default) or a vector of
  two or more values searches a grid and selects by BIC (`n_lambdas`/
  `lambda_min_ratio` control the automatic grid, or a vector of your own
  values is used directly); a single non-negative number fits just that
  lambda, no search. Every other argument (`intercept`, `precision`,
  `penalty_factor`, `standardize_x`, `refit`, `verbose`, ...) is unchanged.
* The returned object is always classed `"EMLMLasso"` and has the same
  shape either way -- `coefficients`, `sigma2`, `D`, `lambda`, `BIC`, etc.
  for the selected (or only) model, plus `search_performed` (`TRUE` if a
  grid was searched), and `lambda_grid`/`bics`/`active_vars`/
  `grid_converged`/`all_models`, which are `NULL` when `search_performed`
  is `FALSE`. 
* `print()`, `summary()`, `coef()`, `fitted()`, `residuals()` and `plot()`
  now dispatch on this one class; `plot()` errors informatively if called
  on a single-lambda fit, since there is no grid to plot.
* `EMLasso_LMM()` and `EMLasso_BIC_Search()` are no longer exported (they
  remain as internal engines called by `EMLMLasso()`); code calling either
  by name needs to switch to `EMLMLasso()`.
* Note: `precision`'s default is now `1e-4` regardless of mode. Previously
  a single-lambda `EMLasso_LMM()` call defaulted to `1e-6` while
  `EMLasso_BIC_Search()` defaulted to `1e-4`; `EMLMLasso()` uses one
  default for both. Pass `precision = 1e-6` explicitly to reproduce the
  old single-lambda default exactly.

## New `intercept` argument; no intercept fit by default (changes default numeric output)

* `EMLasso_LMM()` and `EMLasso_BIC_Search()` gain an `intercept` argument
  (default `FALSE`). It always wins over whatever `formFixed` itself
  implies -- `y ~ x1 + x2` (an implicit intercept) and `y ~ x1 + x2 - 1` (an
  explicit no-intercept) are now treated identically, and `intercept` is
  the single place that controls whether the fixed-effects design matrix
  gets an intercept column.

## `fitted()`/`residuals()` gain a `refit` argument

* `fitted.EMLasso_LMM()` and `residuals.EMLasso_LMM()` (and
  `fitted.EMLasso_BIC_Search()`) gain a `refit` argument (default
  `FALSE`), matching `coef()`'s existing one. `refit = TRUE` returns the
  post-selection `lme4::lmer()` refit's own fitted values (its fixed
  *and* random effects combined coherently by `lme4` itself, rather than
  mixing the refit's fixed effects with the original EM fit's own
  random-effects predictions) or `y` minus those fitted values; the fit
  must have been created with `refit = TRUE`, or these error informatively.

## New `formFixed = y ~ .` shorthand

* `EMLasso_LMM()` and `EMLasso_BIC_Search()` now accept `y ~ .` for
  `formFixed`, meaning "every other column of `data` is a candidate
  fixed-effect covariate", except the response (`y`) and `groupVar`
  itself. A variable used only in `formRandom` (e.g. a random slope) is
  swept in like any other column if it's also wanted as a fixed effect --
  it is not excluded automatically, matching how it would be written out
  by hand. Errors informatively if the expansion leaves no covariates at
  all (e.g. `data` has only the response and `groupVar`).
* The returned/stored `formFixed` (on both `"EMLasso_LMM"` and
  `"EMLasso_BIC_Search"` objects) is now always the *resolved* formula --
  e.g. `y ~ time + x1 + x2 + ...`.

## `print()` now shows only non-zero/selected coefficients

* `print.EMLasso_LMM()`'s "Fixed effects" section now shows only the
  non-zero (lasso-selected) coefficients, matching what `summary()`
  already did -- previously it printed every coefficient, zeros included,
  which for a wide model (e.g. the 101-covariate riboflavin example) made
  the selected model hard to read at a glance.
* When a `refit` is attached, its coefficients are now shown filtered to
  the covariates actually in the refit (dropping the `NA` entries for
  covariates the lasso did not select), again matching `summary()`.
  Previously the full-length, mostly-`NA` vector was printed, which could
  look like every covariate was included in the unpenalized refit when in
  fact only the lasso-selected ones were (`refit_formula` on the fit
  object always reflected the correct, smaller set; only the `print()`
  display was misleading).
* `print.EMLasso_BIC_Search()` inherits this fix automatically (it prints
  `best_model` via `print.EMLasso_LMM()`).

## New post-selection refit (`refit` argument)

* `EMLasso_LMM()` and `EMLasso_BIC_Search()` gain a `refit` argument
  (default `FALSE`). Lasso coefficients are shrunk toward zero by
  construction, which biases them as point estimates even once the right
  covariates are selected. `refit = TRUE` follows the lasso fit with an
  ordinary, unpenalized `lme4::lmer()` fit on just the covariates the
  lasso selected (plus an always-included intercept and the same
  random-effects structure), and attaches it as `refit_fit` (the `merMod`
  object), `refit_coefficients` (`NA` for any covariate the lasso did not
  select), `refit_formula` and `refit_singular`. The original lasso
  `coefficients` are always kept alongside, unaffected.
  `EMLasso_BIC_Search()` refits only the final best-lambda model, not
  every lambda in the grid.
* **The refit is for point estimates only.** Its standard errors,
  p-values and confidence intervals are *not* valid inference: they treat
  the selected covariate set as fixed in advance, when it was in fact
  chosen using the same data (the classic post-selection/selective
  inference problem). This package never computes or reports refit-based
  inferential statistics, and `print()`/`summary()` repeat this caveat
  wherever refit results are shown.
* `coef.EMLasso_LMM()`/`coef.EMLasso_BIC_Search()` gain a `refit`
  argument (default `FALSE`) to retrieve `refit_coefficients` instead of
  the lasso `coefficients`; both error informatively if no refit was
  attached.
* Requires the `lme4` package, only when `refit = TRUE` is actually used.
* The `EMLMLasso-riboflavin` vignette now fits with `refit = TRUE` and
  reports the refit coefficients as its primary results, alongside the
  raw lasso estimates for comparison.

## New real-data example (`riboflavinV100`)

* Added `riboflavinV100`, a longitudinal version of the riboflavin
  production dataset for *B. subtilis* (Buhlmann, Kalisch and Meier, 2014),
  reconstructed from the original microarray file-name metadata to
  reproduce the real-data application in Oliveira, Schumacher and Lachos
  (2026, Section 5): 28 strains, 2-4 observations each, 71 total, with the
  same 100-gene reduction (the genes of largest empirical variance) used
  by the paper itself. See `?riboflavinV100` for full provenance.
* Added a dedicated vignette, `vignette("EMLMLasso-riboflavin")`, fitting
  `EMLasso_BIC_Search()` to `riboflavinV100`.

Initial package release, built around the March 2026 revision of the EM
algorithm (previously distributed as loose scripts: `main-EMLMLasso.R` and
`EMLMwBICsearch-9-11.R`). This is the version to use going forward; the
sections below record how it differs from the earlier scripts, for anyone
comparing results against older output.

## Lambda selection: BIC search replaces cross-validation

* `EMLasso_BIC_Search()` is now the primary way to choose the lasso penalty.
  It fits the EM algorithm over a data-driven grid of lambda values (from an
  initial `glmnet` fit) and selects the value that minimizes BIC.
* This **replaces** the old `EMLasso_cvLMM()`, which re-ran `cv.glmnet` on
  the random-effects-adjusted response at every EM iteration. That approach
  was slow (one k-fold cross-validation per iteration, for every candidate
  model) and, in simulation, produced a lambda path that was not monotone in
  the number of active covariates. `EMLasso_cvLMM()` is not carried
  over to this package; there is no cv-based alternative.

## Core EM algorithm (`EMLasso_LMM`)

* Variance-floor protection: `sigma2` and the M-step's `sigma2` update are
  both clamped at a minimum of `1e-4`, preventing the degenerate
  near-zero-variance solutions that could occur in the earlier script.
* Convergence is now judged on the Euclidean distance between successive
  parameter vectors `(beta, sigma2, D)` (`precision` compared to
  `sqrt(sum((theta[t] - theta[t-1])^2))`), rather than relative change in the
  log-likelihood (`abs(1 - loglik[t-1]/loglik[t])`).
* `penalty_factor` is exposed as an argument, so individual coefficients can
  be exempted from the penalty; by default, an intercept (if present in
  `formFixed`) is left unpenalized and everything else is penalized equally.

## New interface: formulas + data, not raw matrices

* Every function now takes `data`, `formFixed`, `groupVar` and `formRandom`
  arguments (in that order) instead of separate `y`, `x`, `z`, `nj` arguments.
  Rows are reordered internally by `groupVar`, so `data` does not need to
  be pre-sorted by subject.
* `simulate_lasso_lmm()` is a new function that generates example data
  compatible with this interface (used throughout the documentation, the
  vignette and the test suite).

## Documentation and testing

* Added roxygen documentation and runnable `@examples` for every exported
  function, an introductory vignette, and a `testthat` test suite (none of
  which existed for the original scripts).
