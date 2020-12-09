test_that("RawGen", {
  # univariate log distribution
  set.seed(123)
  testthat::expect_true(suppressWarnings(round(RawGen(baboon.parms_df,
  dist = "log")[1, 6][[1]], 2)) == 70.93)
  # univariate truncated distribution
  set.seed(123)
  testthat::expect_true(suppressWarnings(round(RawGen(baboon.parms_df,
  dist = "trunc")[1, 6][[1]], 2)) == 73.62)
  # multivariate distribution
  set.seed(123)
  testthat::expect_true(suppressWarnings(round(RawGen(
  baboon.parms_list)[1, 6][[1]], 2)) == 58.55)
  set.seed(123)
  testthat::expect_true(suppressWarnings(round(RawGen(baboon.parms_df,
  R.res = R)[1, 6][[1]], 2)) == 58.55)
  testthat::expect_error(suppressWarnings(RawGen(baboon.parms_df,
  R.res = R, dist = "log")))
  testthat::expect_error(suppressWarnings(RawGen(baboon.parms_df,
  R.res = R, dist = "log")))
  testthat::expect_error(suppressWarnings(RawGen(R)))
  testthat::expect_error(suppressWarnings(RawGen(baboon.parms_df,
  complete_cases = 95)))
  testthat::expect_error(suppressWarnings(RawGen(baboon.parms_df,
  Pop = 500)))
  testthat::expect_error(suppressWarnings(RawGen(baboon.parms_df,
  Parms = 500)))
  testthat::expect_true(ncol(suppressWarnings(RawGen(baboon.parms_df,
  format = "long"))) == 4)
  testthat::expect_error(suppressWarnings(RawGen(baboon.parms_df,
                                                 R.res = baboon.parms_list)))
  testthat::expect_error(suppressWarnings(RawGen(baboon.parms_list[1:5])))
})
