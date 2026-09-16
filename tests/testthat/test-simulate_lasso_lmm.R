test_that("simulate_lasso_lmm returns the expected shape", {
  dat <- simulate_lasso_lmm(n = 10, nj = 3, p = 4, pstar = 2, seed = 1)
  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 30)
  expect_equal(ncol(dat), 3 + 4) # id, time, y, x1..x4
  expect_true(all(c("id", "time", "y", "x1", "x2", "x3", "x4") %in% names(dat)))
  expect_equal(length(unique(dat$id)), 10)
  expect_equal(as.integer(table(dat$id)), rep(3L, 10))
})

test_that("simulate_lasso_lmm is reproducible with seed and validates pstar", {
  dat1 <- simulate_lasso_lmm(n = 5, nj = 2, p = 3, pstar = 1, seed = 99)
  dat2 <- simulate_lasso_lmm(n = 5, nj = 2, p = 3, pstar = 1, seed = 99)
  expect_equal(dat1, dat2)

  expect_error(simulate_lasso_lmm(n = 5, nj = 2, p = 3, pstar = 5),
               "pstar")
})
