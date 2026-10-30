#' Longitudinal riboflavin production data (100-gene reduction)
#'
#' Gene expression and riboflavin production data on *Bacillus subtilis*,
#' with a longitudinal (repeated-measures) structure recovered from the
#' original microarray metadata, reproducing the real-data application in
#' Oliveira, Schumacher and Lachos (2026), Section 5.
#'
#' @format A data frame with 71 rows (observations) and 103 columns.
#' `strain` (factor, 28 levels) is the *B. subtilis* strain/isolate a
#' given observation belongs to (the grouping/cluster variable), each
#' measured 2 to 4 times; `time` (numeric) is the time of measurement in
#' hours after induction (24, 30, 48 or 52); `y` (numeric) is
#' log-transformed riboflavin (vitamin B2) production rate, the response
#' variable. The remaining 100 columns are numeric log-expression levels,
#' one per gene (column names are the genes' probe identifiers, e.g.
#' \code{YCIC_at}) -- the 100 genes of largest empirical variance across
#' these 71 observations (see Details).
#'
#' @details
#' The expression matrix and response are the standard `riboflavin`
#' dataset (Bühlmann, Kalisch and Meier, 2014), most commonly distributed
#' via the `hdi` package, which records 71 microarray samples and 4088
#' gene expression levels but exposes no grouping or time variable. The
#' original .CEL microarray file names (preserved as that object's row
#' names) encode the *B. subtilis* strain (e.g. `"Fbat1070"`) and the
#' measurement time in hours (e.g. `"PT24"`) each sample belongs to.
#' Parsing these names recovers exactly the longitudinal structure
#' described in Oliveira, Schumacher and Lachos (2026, Section 5): 28
#' strains, each measured 2-4 times, 71 observations in total.
#'
#' Following prior work, genes are filtered to the 100 of largest
#' empirical variance across these 71 observations -- the same
#' `riboflavinV100` reduction used by the paper itself
#' (following Schelldorfer, Bühlmann and van de Geer, 2011), so that the
#' analysis setting matches the paper's (results differ from its Table 5
#' because of the change in \eqn{\lambda} selection; see [EMLMLasso()]) and
#' fitting is fast. See \code{vignette("EMLMLasso-riboflavin")} for a worked
#' example fitting `EMLMLasso()` to it.
#'
#' @references
#' Oliveira DCR, Schumacher FL, Lachos VH (2026). "The use of the EM
#' algorithm for regularization problems in high-dimensional linear
#' mixed-effects models." *Statistical Methods in Medical Research*,
#' 35(2), 240-254. \doi{10.1177/09622802251399913}
#'
#' Bühlmann P, Kalisch M, Meier L (2014). "High-Dimensional Statistics
#' with a View Toward Applications in Biology." *Annual Review of
#' Statistics and Its Application*, 1, 255-278.
#'
#' Schelldorfer J, Bühlmann P, van de Geer S (2011). "Estimation for
#' high-dimensional linear mixed-effects models using
#' \eqn{\ell_1}-penalization." *Scandinavian Journal of Statistics*, 38(2),
#' 197-214.
#'
#' @source
#' Derived from the `riboflavin` dataset distributed with the `hdi`
#' package (\url{https://github.com/cran/hdi}), combined with strain/time
#' metadata parsed from that dataset's original microarray file names.
#'
#' @examples
#' data(riboflavinV100)
#' dim(riboflavinV100)
#' table(table(riboflavinV100$strain))  # group sizes: 2 to 4 per strain
#' # see vignette("EMLMLasso-riboflavin") for a full worked example
"riboflavinV100"
