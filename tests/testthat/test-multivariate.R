test_that("multivariate", {
  testthat::expect_true(multivariate(baboon.parms_list)$p.value[1]==0.2108)
})
