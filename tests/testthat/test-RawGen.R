test_that("RawGen", {
  testthat::expect_true(nrow(RawGen(baboon.parms_df))==604)
})
