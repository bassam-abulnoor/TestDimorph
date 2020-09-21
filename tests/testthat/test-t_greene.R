test_that("t_greene", {
  testthat::expect_true(
    t_greene(
      baboon.parms_df[1:3, ],
      Pop = 2,
      padjust = "none",
      letters = TRUE,
      plot = TRUE,
      es = TRUE,
    )[[1]]$p.value[1] ==
      0.8564
  )
  testthat::expect_error(
    t_greene(
      baboon.parms_df[1, ]
  ))
  testthat::expect_error(
    t_greene(
      baboon.parms_df[1:3, ],Pop = 500
    ))
  testthat::expect_error(
    t_greene(
      baboon.parms_df[1:3, ],sig.level  = 500,Pop = 2
    ))
  testthat::expect_error(
    t_greene(
      R
    ))
  testthat::expect_true(
    t_greene(
      baboon.parms_df[1:3, ],
      Pop = 2,
      padjust = "none",
      alternative = "less"
    )$p.value[1] ==
      0.5718
  )
  testthat::expect_true(
    t_greene(
      baboon.parms_df[1:3, ],
      Pop = 2,
      padjust = "none",
      alternative = "great"
    )$p.value[1] ==
      0.4282
  )
  testthat::expect_true(
    t_greene(baboon.parms_df[1:3,],Pop = 2,es = TRUE)[9][[1]][1]==0.0528)
})
