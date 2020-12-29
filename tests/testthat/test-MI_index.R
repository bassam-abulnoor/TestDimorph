test_that("MI_index", {
  testthat::expect_true(MI_index(Cremains_measurements[1, ])$MI == 0.2008)
  testthat::expect_true(round(MI_index(Cremains_measurements[1, ], rand = F, B = 100)[[4]], 4)
  == 0.2905)
  testthat::expect_true(MI_index(Cremains_measurements[1, ], index_type = "NI")[[2]] == 0.4016)
  expect_doppelganger <- function(title, fig, path = NULL, ...) {
    testthat::skip_if_not_installed("vdiffr")
    vdiffr::expect_doppelganger(title, fig, path = path, ...)
  }
  expect_doppelganger(
    title = "MI_index",
    fig = MI_index(Cremains_measurements[1, ], plot = TRUE,fill="both")[[2]]
  )
})
