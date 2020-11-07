test_that("multivariate", {
  testthat::expect_true(multivariate(baboon.parms_list, es = TRUE)$p.value[1] == 0.2108)
  testthat::expect_true(multivariate(baboon.parms_df, R.res = R, es = TRUE)$p.value[1] == 0.2108)
  testthat::expect_true(multivariate(baboon.parms_list, univariate = TRUE, padjust = "none", es = TRUE)[[2]][[1]]$p.value[1] == 0.6498)
  testthat::expect_error(multivariate(R))

})
