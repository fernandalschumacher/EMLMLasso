# legacy/

These are the original, pre-package development scripts, kept here for
provenance/comparison. They are not part of the installable `EMLMLasso`
package (this folder is excluded via `.Rbuildignore`) and are not sourced
by anything in `R/`.

* `main-EMLMLasso.R` -- original script: `logvero()`, `EM_LMM()` (unpenalized),
  `EMLasso_LMM()` and `EMLasso_cvLMM()` (cv.glmnet-based lambda selection).
* `EMLMwBICsearch-9-11.R` -- the March 2026 revision this package is built
  from: `logveropen()`, `EMLasso_LMM_testv5()`, and the new
  `EMLasso_BIC_Search()`.
* `application.R` -- a data-application script comparing this method against
  `splmm` and `glmmLasso` on simulated data. It sources `main-EMLMLasso.R`
  and will not run as-is against the package; see `NEWS.md` at the package
  root for what changed and `vignette("EMLMLasso-intro")` for the current,
  supported usage.

See `../NEWS.md` for exactly what changed between these scripts and the
package's `EMLasso_LMM()` / `EMLasso_BIC_Search()`.
