test_that("univariate", {
  testthat::expect_true(univariate(baboon.parms_df[1:3, ], Pop = 2)$p.value[1] ==
                          0.6498)
})
