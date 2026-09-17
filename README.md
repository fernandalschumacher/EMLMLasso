# EMLMLasso

An R package for fitting linear mixed models with a lasso penalty on the
fixed effects, for clustered/longitudinal data, via an EM algorithm.

## Installation

```r
# from GitHub
# install.packages("remotes")
remotes::install_github("fernandalschumacher/EMLMLasso", build_vignettes = TRUE)
```

## Usage

```r
library(EMLMLasso)

dat <- simulate_lasso_lmm(n = 40, nj = 5, p = 10, pstar = 3, seed = 2026)
fixedForm <- as.formula(paste("y ~", paste0("x", 1:10, collapse = " + ")))

# select lambda by BIC over a data-driven grid (the default, and the
# recommended way to call this for most uses)
search <- EMLMLasso(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time)
search
coef(search)
plot(search)

# or fit a single lambda directly, by passing one in
fit <- EMLMLasso(dat, formFixed = fixedForm, groupVar = "id", formRandom = ~time,
                  lambda = 0.15)
```

See `vignette("EMLMLasso-intro")` for a full walkthrough, and `NEWS.md` for
what changed from the earlier development scripts.

## Citation

This package implements the method introduced in:

> Oliveira DCR, Schumacher FL, Lachos VH (2026). The use of the EM algorithm
> for regularization problems in high-dimensional linear mixed-effects
> models. *Statistical Methods in Medical Research*, 35(2), 240-254.
> [doi:10.1177/09622802251399913](https://doi.org/10.1177/09622802251399913)

Run `citation("EMLMLasso")` for a BibTeX-ready entry.

## Acknowledgments

Portions of this package's documentation were drafted with the assistance
of Claude (Anthropic).
