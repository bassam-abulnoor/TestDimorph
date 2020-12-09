library(TestDimorph)
if (requireNamespace('testthat', quietly = TRUE)) {
  library('testthat')
} else {
  message("'testthat' not available")
}
test_check("TestDimorph")
