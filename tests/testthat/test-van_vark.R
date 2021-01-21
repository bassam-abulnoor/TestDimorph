test_that("van_vark", {
  library(TestDimorph)
  testthat::expect_error(TestDimorph::van_vark(Howells_summary, W = Howells_V, Trait = 500, Pop = 500))
  testthat::expect_error(TestDimorph::van_vark(Howells_summary, W = Howells_V, plot = 55))
  testthat::expect_warning(TestDimorph::van_vark(Howells_summary, W = Howells_V, plot = TRUE, q = 1))
  testthat::expect_true(TestDimorph::van_vark(Howells_summary, W = Howells_V)[[1]][[4]][[1]] == 0.5271)
})
