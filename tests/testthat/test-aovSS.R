test_that("aovSS", {
  testthat::expect_true(suppressWarnings(round(
    aovSS(
      baboon.parms_df[1:3, ],
      Pop = 2,
      digits = 3,
      es = TRUE,
      letters = TRUE,
      pairwise = TRUE,
    )[[1]]$p.value[1],
    3
  )) == 0.325)
})
