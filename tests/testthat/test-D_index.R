test_that("D_index", {
  library(TestDimorph)
  testthat::expect_true(round(D_index(Cremains_measurements[1, ])$D, 4) == 0.5983)
  testthat::expect_true(round(D_index(Cremains_measurements[1, ], rand = F, B = 100)[[4]], 4)
  == 0.7701)
<<<<<<< HEAD
  expect_doppelganger <- function(title, fig, path = NULL, ...) {
    testthat::skip_if_not_installed("vdiffr")
    vdiffr::expect_doppelganger(title, fig, path = path, ...)
  }
  expect_doppelganger(
    title = "D_index",
    fig = D_index(Cremains_measurements[1, ], plot = TRUE, fill = "both")[[2]]
  )
=======
>>>>>>> new5
})
