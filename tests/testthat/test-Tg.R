test_that("Tg", {
  testthat::expect_true(
    suppressWarnings(Tg(
      baboon.parms_df[1:3, ],
      Pop = 2,
      padjust = "none",
      letters = TRUE,
      plot = TRUE,
      es = TRUE,
    )[[1]]$p.value[1] ==
      0.8564)
  )
  testthat::expect_error(
    suppressWarnings(Tg(
      baboon.parms_df[1, ]
    )))
  testthat::expect_error(
    suppressWarnings(Tg(
      baboon.parms_df[1:3, ],Pop = 500
    )))
  testthat::expect_error(
    suppressWarnings(Tg(
      baboon.parms_df[1:3, ],sig.level  = 500,Pop = 2
    )))
  testthat::expect_error(
    suppressWarnings(Tg(
      R
    )))
  testthat::expect_true(
    suppressWarnings(Tg(
      baboon.parms_df[1:3, ],
      Pop = 2,
      padjust = "none",
      alternative = "less"
    )$p.value[1] ==
      0.5718
  ))
  testthat::expect_true(
    suppressWarnings(Tg(
      baboon.parms_df[1:3, ],
      Pop = 2,
      padjust = "none",
      alternative = "great"
    )$p.value[1] ==
      0.4282
  ))
})
