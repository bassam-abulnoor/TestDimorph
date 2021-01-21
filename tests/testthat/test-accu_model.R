test_that("accu_model", {
  testthat::skip_if_not_installed("e1071")
  testthat::expect_true(round(
    accu_model(
      f = Sex ~ GOL + NOL + BNL,
      x = Howells,
      y = Howells,
      byPop = FALSE
    )[[2]][[3]][[1]],
    3
  ) == 0.789)

  testthat::expect_true(round(
    accu_model(
      f = Sex ~ GOL + NOL + BNL,
      x = Howells,
      y = Howells,
      byPop = TRUE,
      Pop = 2
    )[[2]][[1]][[3]][[1]],
    3
  ) == 0.91)
  set.seed(123)
  testthat::expect_true(round(
    accu_model(
      f = Sex ~ GOL + NOL + BNL,
      x = Howells,
      byPop = FALSE
    )[[2]][[3]][[1]],
    3
  ) == 0.811)
  set.seed(123)
  testthat::expect_true(round(
    accu_model(
      f = Sex ~ GOL + NOL + BNL,
      x = Howells,
      byPop = TRUE,
      Pop = 2
    )[[2]][[3]][[3]][[1]],
    3
  ) == 0.762)
  testthat::expect_error(
    accu_model(
      f = Sex ~ GOL + NOL + BNL,
      x = matrix(NA),
      byPop = 50,
      Pop = 200,
      y = matrix(NA),
      plot = 98,
      Sex = 500,
      post. = "ll",
      ref. = "kl"
    )
  )
  testthat::expect_warning(accu_model(
    f = Sex ~ GOL + NOL + BNL,
    x = Howells,
    byPop = TRUE,
    Pop = NULL
  ))
})
