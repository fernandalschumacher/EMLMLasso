#' EMLMLasso: EM Algorithm for Lasso-Penalized Linear Mixed Models
#'
#' Fits linear mixed models with a lasso penalty on the fixed effects, for
#' clustered/longitudinal data, via an EM algorithm. See [EMLMLasso()], the
#' package's single entry point: pass a single `lambda` to fit one penalty
#' value directly, or leave it `NULL` (or supply a vector) to search a grid
#' and select by BIC. [simulate_lasso_lmm()] generates example data with a
#' sparse true fixed-effects vector.
#'
#' @references
#' Oliveira DCR, Schumacher FL, Lachos VH (2026). The use of the EM algorithm
#' for regularization problems in high-dimensional linear mixed-effects
#' models. *Statistical Methods in Medical Research*, 35(2), 240-254.
#' \doi{10.1177/09622802251399913}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
