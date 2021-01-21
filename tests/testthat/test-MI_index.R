test_that("MI_index", {
  library(TestDimorph)
  testthat::expect_true(MI_index(Cremains_measurements[1, ])$MI == 0.2008)
  testthat::expect_true(round(MI_index(Cremains_measurements[1, ], rand = F, B = 100)[[4]], 4)
  == 0.2905)
  testthat::expect_true(MI_index(Cremains_measurements[1, ], index_type = "NI")[[2]] == 0.4016)
})
