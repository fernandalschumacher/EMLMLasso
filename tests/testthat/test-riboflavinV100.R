test_that("riboflavinV100 has the expected longitudinal structure", {
  data(riboflavinV100, package = "EMLMLasso", envir = environment())

  expect_s3_class(riboflavinV100, "data.frame")
  expect_equal(nrow(riboflavinV100), 71)
  expect_true(is.factor(riboflavinV100$strain))
  expect_equal(nlevels(riboflavinV100$strain), 28)

  grp_sizes <- table(riboflavinV100$strain)
  expect_equal(sum(grp_sizes), 71)
  expect_true(all(grp_sizes %in% 2:4))

  expect_true(all(c("time", "y") %in% names(riboflavinV100)))
  expect_false(anyNA(riboflavinV100))

  gene_cols <- setdiff(names(riboflavinV100), c("strain", "time", "y"))
  expect_equal(length(gene_cols), 100)
  expect_true(all(vapply(riboflavinV100[gene_cols], is.numeric, logical(1))))
})
