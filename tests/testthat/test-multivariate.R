test_that("multivariate", {
  library(TestDimorph)
  testthat::expect_true(multivariate(baboon.parms_list)$p.value[3] == 0.2108)
  testthat::expect_true(multivariate(baboon.parms_df, R.res = baboon.parms_R)$p.value[3] ==
    0.2108)
  testthat::expect_true(
    multivariate(
      baboon.parms_list,
      univariate = TRUE,
      padjust = "none"
    )[[2]][[1]]$p.value[3] == 0.6498
  )
  testthat::expect_error(multivariate(baboon.parms_R))
  testthat::expect_true(multivariate(baboon.parms_list, es_manova = "eta")
  [[9]][1] == 0.0656)
  testthat::expect_true(length(multivariate(
    baboon.parms_df,
    R.res = baboon.parms_R,
    univariate = T
  )[[2]]) == 4)
  testthat::expect_error(multivariate(baboon.parms_list, univariate = 55))
  testthat::expect_error(multivariate(Howells))
  testthat::expect_error(multivariate(baboon.parms_df,
    R.res = baboon.parms_R,
    Pop = 800
  ))
  testthat::expect_error(multivariate(baboon.parms_df,
    R.res = baboon.parms_R,
    Trait = 800
  ))
  testthat::expect_error(multivariate(baboon.parms_df))
  testthat::expect_error(multivariate(baboon.parms_df, R.res = Howells))
  testthat::expect_error(multivariate(list("A", "B")))
  testthat::expect_true(multivariate(Howells_summary, Howells_R,
    manova_test_statistic =
      "P"
  )$p.value[3] == 0.0709)
  testthat::expect_true(multivariate(Howells_summary, Howells_R,
    manova_test_statistic =
      "HL"
  )$p.value[3] == 0.0681)
  testthat::expect_true(multivariate(Howells_summary, Howells_R,
    manova_test_statistic =
      "R"
  )$p.value[3] == 0.0032)
  testthat::expect_true(multivariate(Howells_summary, Howells_R,
    manova_test_statistic =
      "R", es_manova = "eta"
  )$eta[3] == 0.0522)
  testthat::expect_true(multivariate(Howells_summary, Howells_R,
    manova_test_statistic =
      "P", es_manova = "eta"
  )$eta[3] == 0.0265)
  testthat::expect_true(multivariate(Howells_summary, Howells_R,
    manova_test_statistic =
      "HL", es_manova = "eta"
  )$eta[3] == 0.0269)
  testthat::expect_true(multivariate(baboon.parms_df, baboon.parms_R,
    univariate = T
  )[[2]][[1]]$p.value[3] == 0.6498)
  testthat::expect_true(multivariate(Howells_summary, Howells_R,
    manova_test_statistic = "HL", es_manova = "eta", interact_manova = FALSE,
    type_manova = "I"
  )[[1]]$p.value[1] == 0)
  testthat::expect_true(multivariate(Howells_summary, Howells_R,
    manova_test_statistic = "R", es_manova = "eta", interact_manova = FALSE,
    type_manova = "II"
  )$p.value[1] == 0)
  testthat::expect_true(multivariate(Howells_summary, Howells_R,
    manova_test_statistic = "P",
    es_manova = "eta", interact_manova = TRUE, type_manova = "I"
  )[[1]]$p.value[3] == 0.0709)
  testthat::expect_true(multivariate(Howells_summary, Howells_R,
    manova_test_statistic = "HL",
    es_manova = "eta", interact_manova = TRUE, type_manova = "II"
  )$p.value[3] == 0.0681)
  testthat::expect_true(multivariate(Howells_summary, Howells_R,
    manova_test_statistic = "W",
    es_manova = "eta", interact_manova = TRUE, type_manova = "III"
  )$p.value[3] == 0.0695)
  testthat::expect_error(multivariate(baboon.parms_list,
    type_manova = "III",
    interact_manova = FALSE
  ))
  testthat::expect_true(length(multivariate(baboon.parms_list, type_manova = "I")) == 2)
})
