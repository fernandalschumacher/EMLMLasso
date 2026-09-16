## Build the longitudinal riboflavin example dataset bundled with this package.
##
## Source: the log-expression matrix and response are the standard
## `riboflavin` dataset (Buhlmann, Kalisch & Meier, 2014, "High-dimensional
## statistics with a view toward applications in biology", Annual Review of
## Statistics and Its Application, distributed e.g. via the `hdi` package,
## https://github.com/cran/hdi). That object has no explicit grouping
## variable, but its row names are the original microarray (.CEL) file
## names, which encode the B. subtilis strain ("Fbat<id>") and the
## measurement time in hours ("PT<hours>") -- e.g.
## "knh_490_BS0001_Fbat1070_PT24_1.Rep.CEL" is strain Fbat1070 at 24 hours.
## Parsing these recovers exactly the longitudinal structure described in
## Oliveira, Schumacher & Lachos (2026, Section 5): 28 strains measured
## 2-4 times over 96 hours, 71 observations total (n_i in {2,3,4}).
##
## Genes are filtered to the top 100 by empirical variance across these 71
## observations, matching the paper's own `riboflavinV100` reduction (also
## used by `hdi`, `DLFM`, `qut`, `FPCdpca`, etc.) -- the standard,
## well-known reduction of this dataset, so results are directly comparable
## to the paper's own Table 5 and the vignette's BIC search runs quickly.
stopifnot(requireNamespace("utils", quietly = TRUE))

## `riboflavin_raw.csv` has columns y, x.<GENE_at> (4089 cols) with rows in
## the same order as hdi::riboflavin (verified against a live session).
## It carries no row names, so the strain/time metadata below is supplied
## separately, in that same row order.
raw <- read.csv("data-raw/riboflavin_raw.csv", check.names = FALSE)

cel_names <- c(
  "b_Fbat107PT24.CEL", "b_Fbat107PT30.CEL", "b_Fbat107PT48.CEL", "b_Fbat107PT52.CEL",
  "knh_102_Fbat289PT24.CEL", "knh_103_Fbat289PT48.CEL",
  "knh_490_BS0001_Fbat1070_PT24_1.Rep.CEL", "knh_491_BS0001_Fbat1070_PT30_1.Rep.CEL",
  "knh_492_BS0001_Fbat1070_PT48_1.Rep.CEL", "knh_494_BS0001_Fbat1077_PT24_2.Rep.CEL",
  "knh_495_BS0001_Fbat1077_PT30_2.Rep.CEL", "knh_496_BS0001_Fbat1077_PT48_2.Rep.CEL",
  "knh_498_BS0001_Fbat1091_PT24_3.Rep.CEL", "knh_499_BS0001_Fbat1091_PT30_3.Rep.CEL",
  "knh_500_BS0001_Fbat1091_PT48_3.Rep.CEL", "knh_502_BS3416_E_Fbat1071_PT24_1.Rep.CEL",
  "knh_503_BS3416_E_Fbat1071_PT30_1.Rep.CEL", "knh_504_BS3416_E_Fbat1071_PT48_1.Rep.CEL",
  "knh_506_BS3416_E_Fbat1078_PT24_2.Rep.CEL", "knh_507_BS3416_E_Fbat1078_PT30_2.Rep.CEL",
  "knh_508_BS3416_E_Fbat1078_PT48_2.Rep.CEL", "knh_510_BS3416_E_Fbat1092_PT24_3.Rep.CEL",
  "knh_512_BS3416_E_Fbat1092_PT48_3.Rep.CEL", "knh_532_BS3416_E_Fbat1092_PT30_3.Rep_wh.CEL",
  "knh_655_BS5009_Fbat1374_PT24_1.Repl._No.2_repeated RNA.CEL",
  "knh_657_BS5009_Fbat1379_PT24_2.Repl._No.6_repeated RNA.CEL",
  "knh_659_BS5009_Fbat1383_PT24_3.Repl._No.10_repeated RNA.CEL",
  "knh_660_BS5009_Fbat1374_PT30_1.Repl._No.3_rep_cDNA_new Yeast.CEL",
  "knh_661_BS5009_Fbat1374_PT48_1.Repl._No.4_rep_cDNA_new Yeast.CEL",
  "knh_662_BS5009_Fbat1379_PT30_2.Repl._No.7_rep_cDNA_new Yeast.CEL",
  "knh_663_BS5009_Fbat1379_PT48_2.Repl._No.8_rep_cDNA_new Yeast.CEL",
  "knh_664_BS5009_Fbat1383_PT30_3.Repl._No.11_rep_cDNA_new Yeas.CEL",
  "knh_665_BS5009_Fbat1383_PT48_3.Repl._No.12_rep_cDNA_new Yeas.CEL",
  "knh_679_BS5254_Fbat1442_PT24_1.Repl._No.38.CEL", "knh_680_BS5254_Fbat1442_PT30_1.Repl._No.39.CEL",
  "knh_681_BS5254_Fbat1442_PT48_1.Repl._No.40.CEL", "knh_683_BS5254_Fbat1443_PT24_2.Repl._No.42.CEL",
  "knh_684_BS5254_Fbat1443_PT30_2.Repl._No.43.CEL", "knh_685_BS5254_Fbat1443_PT48_2.Repl._No.44.CEL",
  "knh_687_BS5254_Fbat1444_PT24_3.Repl._No.46.CEL", "knh_688_BS5254_Fbat1444_PT30_3.Repl._No.47.CEL",
  "knh_689_BS5254_Fbat1444_PT48_3.Repl._No.48.CEL", "knhb_039_Fbat152PT24.CEL",
  "knhb_040_Fbat152PT48.CEL", "knhb_041_Fbat152PT52.CEL", "knhb_046_Fbat202PT24.CEL",
  "knhb_047_Fbat202PT48.CEL", "knhb_049_Fbat203PT24.CEL", "knhb_050_Fbat203PT48.CEL",
  "knhb_090_Fbat284PT24.CEL", "knhb_091_Fbat284PT48.CEL", "knhb_098_Fbat287PT24.CEL",
  "knhb_099_Fbat287PT48.CEL", "knhb_169_Fbat395PT24.CEL", "knhb_170_Fbat395PT48.CEL",
  "knhb_174_Fbat396PT24.CEL", "knhb_175_Fbat396PT48.CEL", "knhb_180_Fbat397PT24.CEL",
  "knhb_181_Fbat397PT48.CEL", "knhb_186_Fbat398PT24.CEL", "knhb_187_Fbat398PT48.CEL",
  "knhb_237_Fbat394PT24.CEL", "knhb_238_Fbat394PT48.CEL", "knhb_247_Fbat525PT24.CEL",
  "knhb_248_Fbat525PT48.CEL", "knhb_251_Fbat526PT24.CEL", "knhb_252_Fbat526PT48.CEL",
  "knhb_255_Fbat527PT24.CEL", "knhb_256_Fbat527PT48.CEL", "knhb_259_Fbat528PT24.CEL",
  "knhb_260_Fbat528PT48.CEL"
)
stopifnot(length(cel_names) == nrow(raw))

m <- regmatches(cel_names, regexec("Fbat(\\d+)_?PT(\\d+)", cel_names))
strain <- vapply(m, `[`, character(1), 2)
time_h <- as.numeric(vapply(m, `[`, character(1), 3))
stopifnot(!anyNA(time_h), all(nchar(strain) > 0))

## sanity check against Oliveira, Schumacher & Lachos (2026), Section 5
grp_sizes <- table(strain)
stopifnot(length(grp_sizes) == 28, sum(grp_sizes) == 71,
          all(grp_sizes %in% 2:4))

y <- raw[[1]]
x <- as.matrix(raw[, -1, drop = FALSE])
colnames(x) <- make.names(sub("^x\\.", "", colnames(x)), unique = TRUE)

## keep the 100 genes of largest empirical variance across these 71 obs
## (the paper's own riboflavinV100 reduction)
gene_var <- apply(x, 2, stats::var)
keep <- order(gene_var, decreasing = TRUE)[1:100]
x <- x[, sort(keep), drop = FALSE]

riboflavinV100 <- data.frame(
  strain = factor(strain, levels = unique(strain)),
  time = time_h,
  y = y,
  x,
  check.names = FALSE,
  row.names = NULL
)
riboflavinV100 <- riboflavinV100[order(match(riboflavinV100$strain, unique(strain)), riboflavinV100$time), ]
rownames(riboflavinV100) <- NULL

save(riboflavinV100, file = "data/riboflavinV100.rda", compress = "xz")
