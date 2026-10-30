## Resubmission

This is a resubmission. In this version I have:

* Removed the single quotes around function names in the Description
  field of DESCRIPTION (e.g. EMLMLasso()).
* Added \value sections to coef.EMLMLasso.Rd and plot.EMLMLasso.Rd,
  describing the class and meaning of the output (plot() is documented as
  invisibly returning its input, called for its side effect).

In addition, I made the following changes since the previous submission:

* Added a `center` argument to EMLMLasso() (default TRUE), so that the
  response and covariates are centred before fitting when no intercept is
  included.
* Fixed a bug where fits that reached the maximum number of EM iterations
  were reported as converged.
* Clarified in the documentation how this implementation differs from the
  published method, and updated the vignettes accordingly.

## R CMD check results

Local check (Windows 11, R 4.6.0, --as-cran): 0 errors | 0 warnings | 0 notes

On CRAN's incoming checks I expect one NOTE:

* This is a new submission.
