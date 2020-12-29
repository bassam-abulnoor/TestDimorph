test_that("extract_sum", {
  testthat::expect_true(TestDimorph::extract_sum(Howells, test = 4)$p.value[3] == 0.0695)
  testthat::expect_true(TestDimorph::extract_sum(Howells,
    test = 1,
    padjust = "none"
  )$p.value[1] == 0.2708)
  testthat::expect_true(TestDimorph::extract_sum(Howells,
    test = 2,
    padjust = "none"
  )$p.value[3] ==
    0.7243)
  testthat::expect_true(TestDimorph::extract_sum(Howells, test = 3)[[1]][[6]][1] == 0)
  testthat::expect_true(TestDimorph::extract_sum(Howells,
    test = 5,
    plot = F
  )[[4]][[1]] == 0.5271)
})
