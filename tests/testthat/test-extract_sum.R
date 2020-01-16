test_that("extract_sum", {
  testthat::expect_true(TestDimorph::extract_sum(Howells,test = 4)$p.value[1]==0.0695)
})
