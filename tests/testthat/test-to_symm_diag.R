test_that("to_symm_diag", {
  testthat::expect_true(length(to_symm_diag(c(3089, 1079, 785, 574, 351, 350)))==
                          9)
  testthat::expect_true(is.matrix(to_symm_diag(c(3089, 1079, 785, 574, 351, 350))))
})
