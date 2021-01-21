test_that("Hedges_gx", {
  testthat::expect_true(round(Hedges_g(Cremains_measurements[1, ])$g, 4) == 1.6272)
  testthat::expect_true(round(Hedges_g(Cremains_measurements[1, ], rand = F, B = 100)[[2]], 4)
  == 0.8933)
})
