test_that("pMatrix", {
  testthat::expect_true(suppressWarnings(pMatrix(extract_sum(
    Howells,
    test = 1, run = FALSE
  ), padjust = "none")[1, 4] == 0.4388))
  testthat::expect_true(suppressWarnings(pMatrix(extract_sum(
    Howells,
    test = 1, run = FALSE
  ))[1, 4] == 1))
  testthat::expect_true(suppressWarnings(pMatrix(extract_sum(
    Howells,
    test = 1, run = FALSE
  ), padjust = "none", alternative = "less")[1, 4] == 0.7806))

  testthat::expect_true(suppressWarnings(pMatrix(extract_sum(
    Howells,
    test = 1, run = FALSE
  ), padjust = "none", alternative = "great")[1, 4] == 0.2194))

  testthat::expect_error(suppressWarnings(pMatrix(extract_sum(
    Howells,
    test = 1, run = FALSE
  )[1:3])))
  testthat::expect_error(suppressWarnings(pMatrix(R)))
})
