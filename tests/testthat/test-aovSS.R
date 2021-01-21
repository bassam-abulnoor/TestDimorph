test_that("aovSS", {
  library(TestDimorph)
  testthat::expect_true(round(
    aovSS(
      baboon.parms_df[1:3, ],
      Pop = 2,
      digits = 3,
      letters = TRUE,
      pairwise = TRUE,
    )[[1]]$p.value[1],
    3
  ) == 0.325)
  testthat::expect_true(aovSS(
    baboon.parms_df[1:3, ],
    Pop = 2,
    digits = 3,
    letters = TRUE,
    pairwise = TRUE,
    es_anova = "f"
  )$`Female model`[[8]][1] == 0.028)
  testthat::expect_true(aovSS(
    baboon.parms_df[1:3, ],
    Pop = 2,
    digits = 3,
    letters = TRUE,
    pairwise = TRUE,
    es_anova = "eta"
  )$`Female model`[[8]][1] == 0.027)
  testthat::expect_error(aovSS(
    baboon.parms_df[1:3, ],
    Pop = 2,
    digits = 3,
    letters = TRUE,
    pairwise = TRUE,
    es_anova = "qq"
  )$`Female model`[[8]][1] == 0.028)
  testthat::expect_error(
    aovSS(
      x = matrix(NA),
      Pop = 500,
      digits = 3,
      es = 900,
      letters = 900,
      pairwise = 900,
      sig.level = 900
    )
  )
})
