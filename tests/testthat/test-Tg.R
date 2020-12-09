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
    ))
  )
  testthat::expect_error(
    suppressWarnings(Tg(
      baboon.parms_df[1:3, ],
      Pop = 500
    ))
  )
  testthat::expect_error(
    suppressWarnings(Tg(
      baboon.parms_df[1:3, ],
      sig.level = 500, Pop = 2
    ))
  )
  testthat::expect_error(
    suppressWarnings(Tg(
      R
    ))
  )
  testthat::expect_true(
    suppressWarnings(Tg(
      baboon.parms_df[1:3, ],
      Pop = 2,
      padjust = "none",
      alternative = "less"
    )$p.value[1] ==
      0.5718)
  )
  testthat::expect_true(
    suppressWarnings(Tg(
      baboon.parms_df[1:3, ],
      Pop = 2,
      padjust = "none",
      alternative = "great"
    )$p.value[1] ==
      0.4282)
  )
  testthat::expect_true(
    suppressWarnings(Tg(data.frame(
      Pop = c("Ireland", "Colombia"),
      m = c(347, 317),
      M.mu = c(172.9, 163.3),
      M.sdev = c(6.34, 6.11),
      f = c(261, 317),
      F.mu = c(159, 151.3),
      F.sdev = c(5.35, 5.31)
    )))[[7]] == 0.0044
  )
})
