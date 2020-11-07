test_that("AccuModel", {
  testthat::expect_true(suppressWarnings(round(
    AccuModel(
      f = Sex ~ GOL + NOL + BNL,
      x = Howells,
      y = Howells,
      byPop = FALSE
    )[[3]][[1]],
    3
  )) == 0.789)

  testthat::expect_true(suppressWarnings(round(
    AccuModel(
      f = Sex ~ GOL + NOL + BNL,
      x = Howells,
      y = Howells,
      byPop = TRUE,
      Pop = 2
    )[[1]][[3]][[1]],
    3
  )) == 0.91)

  vdiffr::expect_doppelganger(
    title = "Single plot 2",
    fig = suppressWarnings(AccuModel(
      f = Sex ~ GOL + NOL + BNL,
      x = Howells,
      y = Howells,
      byPop = FALSE,
      plot = TRUE
    ))[[1]]
  )
  vdiffr::expect_doppelganger(
    title = "Pop plot 2",
    fig = suppressWarnings(AccuModel(
      f = Sex ~ GOL + NOL + BNL,
      x = Howells,
      y = Howells,
      byPop = TRUE,
      Pop = 2,
      plot = TRUE
    ))[[1]]
  )
  set.seed(123)
  testthat::expect_true(suppressWarnings(round(
    AccuModel(
      f = Sex ~ GOL + NOL + BNL,
      x = Howells,
      byPop = FALSE
    )[[3]][[1]],
    3
  )) == 0.811)
  set.seed(123)
  testthat::expect_true(suppressWarnings(round(
    AccuModel(
      f = Sex ~ GOL + NOL + BNL,
      x = Howells,
      byPop = TRUE,
      Pop = 2
    )[[1]][[3]][[1]],
    3
  )) == 0.909)
})
