test_that("tGreene", {
  testthat::expect_true(Tg(baboon.parms_df[1:2,],Pop = 2)$p.value==0.3413)
})
