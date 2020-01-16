test_that("AccuModel", {
  testthat::expect_true(round(AccuModel(f = Sex~GOL+NOL+BNL,x = Howells,y = Howells,byPop = FALSE)[[3]][[1]],3)==0.789)
})
