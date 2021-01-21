test_that("univariate", {
  testthat::expect_true(univariate(baboon.parms_df[1:3, ], Pop = 2)$p.value[3] ==
    0.6498)
  testthat::expect_true(univariate(baboon.parms_df[1:3, ],
    Pop = 2,
    pairwise = T
  )[[2]][[7]][1] == 0.3413)
  testthat::expect_error(univariate(baboon.parms_df[1:3, ], Pop = 200))
  testthat::expect_error(univariate(baboon.parms_df[1:3, 1], Pop = 200))
  testthat::expect_error(univariate(baboon.parms_df[1, 1], Pop = 200))
  testthat::expect_error(univariate(baboon.parms_df[1, 3], es = 500))
  testthat::expect_error(univariate(baboon.parms_df[1, 3], pairwise = 500))
  testthat::expect_true(univariate(baboon.parms_df[1:10, ],
    es_anova = "f"
  )[[8]][3] == 0.0072)
  testthat::expect_true(univariate(baboon.parms_df[1:10, ], es_anova = "eta")
  [[8]][3] == 0.0071)
  testthat::expect_error(univariate(baboon.parms_df[1:10, ], es_anova = "omega")
  [[8]][3] == 0.0029)
  testthat::expect_error(univariate(baboon.parms_df[1:10, ], es_anova = "qq")
  [[8]][1] == 0.0029)
  testthat::expect_true(
    univariate(data.frame(
      Pop = c("Ireland", "Colombia"),
      m = c(347, 317),
      M.mu = c(172.9, 163.3),
      M.sdev = c(6.34, 6.11),
      f = c(261, 317),
      F.mu = c(159, 151.3),
      F.sdev = c(5.35, 5.31)
    ))[[6]][3] == 0.0044
  )
  testthat::expect_true(length(univariate(baboon.parms_df,
    interact_anova = F, type_anova = "I", es_anova =
      "eta", pairwise = F
  )) == 2)
  testthat::expect_true(length(univariate(baboon.parms_df,
    interact_anova = T, type_anova = "I", es_anova =
      "eta", pairwise = F
  )) == 2)
  testthat::expect_true(univariate(baboon.parms_df,
    interact_anova = F, type_anova = "II", es_anova =
      "eta", pairwise = F
  )$p.value[1] == 0)
  testthat::expect_true(univariate(baboon.parms_df,
    interact_anova = T, type_anova = "II", es_anova =
      "eta", pairwise = F
  )$p.value[1] == 0)
  testthat::expect_error(univariate(baboon.parms_df,
    interact_anova = F, type_anova = "III", es_anova =
      "eta", pairwise = F
  )$p.value[1] == 0)
  testthat::expect_true(univariate(baboon.parms_df,
    interact_anova = T, type_anova = "III", es_anova =
      "eta", pairwise = F
  )$p.value[3] == 0.0755)
})
