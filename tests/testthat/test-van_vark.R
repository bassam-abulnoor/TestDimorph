test_that("van_vark", {
  library(TestDimorph)
  testthat::expect_error(TestDimorph::van_vark(Howells_summary, W = Howells_V, Trait = 500, Pop = 500))
  testthat::expect_error(TestDimorph::van_vark(Howells_summary, W = Howells_V, plot = 55))
  testthat::expect_warning(TestDimorph::van_vark(Howells_summary, W = Howells_V, plot = TRUE, q = 1))
  testthat::expect_true(TestDimorph::van_vark(Howells_summary, W = Howells_V)[[1]][[4]][[1]] == 0.5271)
  expect_doppelganger <- function(title, fig, path = NULL, ...) {
    testthat::skip_if_not_installed("vdiffr")
    vdiffr::expect_doppelganger(title, fig, path = path, ...)
  }
  expect_doppelganger(
    title = "TestDimorph::van_vark",
    fig = TestDimorph::van_vark(Howells_summary, Howells_V, plot = TRUE)[[2]]
  )
})
