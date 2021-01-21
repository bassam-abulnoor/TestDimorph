test_that("t_greene", {
  library(TestDimorph)
  testthat::expect_true(
    t_greene(
      baboon.parms_df[1:3, ],
      Pop = 2,
      padjust = "none",
      letters = TRUE,
      plot = TRUE
    )[[1]]$p.value[1] ==
      0.3413
  )
  testthat::expect_error(
    t_greene(
      baboon.parms_df[1, ]
    )
  )
  testthat::expect_error(
    t_greene(
      baboon.parms_df[1:3, ],
      Pop = 500
    )
  )
  testthat::expect_error(
    t_greene(
      baboon.parms_df[1:3, ],
      CI = 500, Pop = 2
    )
  )
  testthat::expect_error(
    t_greene(
      R
    )
  )
  testthat::expect_true(
    t_greene(
      baboon.parms_df[1:3, ],
      Pop = 2,
      padjust = "none",
      alternative = "less"
    )$p.value[1] ==
      0.8293
  )
  testthat::expect_true(
    t_greene(
      baboon.parms_df[1:3, ],
      Pop = 2,
      padjust = "none",
      alternative = "great"
    )$p.value[1] ==
      0.1707
  )
  testthat::expect_true(
    t_greene(baboon.parms_df[1:3, ], Pop = 2, es = "d")[9][[1]][1] == 0.0824
  )
  testthat::expect_true(
    t_greene(baboon.parms_df[1:3, ], Pop = 2, es = "g")[9][[1]][1] == 0.0823
  )
  testthat::expect_error(
    t_greene(baboon.parms_df[1:3, ], Pop = 2, es = "qq")[9][[1]][1] == 0.0181
  )
  testthat::expect_true(
    t_greene(data.frame(
      Pop = c("Ireland", "Colombia"),
      m = c(347, 317),
      M.mu = c(172.9, 163.3),
      M.sdev = c(6.34, 6.11),
      f = c(261, 317),
      F.mu = c(159, 151.3),
      F.sdev = c(5.35, 5.31)
    ))[[7]] == 0.0044
  )
  testthat::expect_error(
    t_greene(
      baboon.parms_list
    )
  )
  testthat::expect_error(
    t_greene(
      Howells
    )
  )
  testthat::expect_error(
    t_greene(
      baboon.parms_df[1:3, ]
    )
  )
  testthat::expect_error(
    suppressWarnings(t_greene(
      baboon.parms_df[1:3, ],
      letters = 55
    ))
  )
  testthat::expect_error(
    suppressWarnings(t_greene(
      baboon.parms_df[1:3, ],
      plot = 55
    ))
  )
})
