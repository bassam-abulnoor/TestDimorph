test_that("D_index", {
  library(TestDimorph)
  testthat::expect_true(round(D_index(Cremains_measurements[1, ])$D, 4) == 0.5983)
  testthat::expect_true(round(D_index(Cremains_measurements[1, ], rand = F, B = 100)[[4]], 4)
  == 0.7701)
})
