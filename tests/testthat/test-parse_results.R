context("parse_results")

test_that("Identify blocks of large covariance matrices", {
  
  ## airemlf90_log_4.txt - 10-variate model.
  ## with large residual covariance matrix (10x10) which wraps lines
  reml_log <- readLines(file.path(testdata, "airemlf90_log_4.txt"))
  resid_var_ini_line <- grep("Residual variance", reml_log)
  
  expect_error(
    bl <- extract_block(resid_var_ini_line + 1, reml_log),
    NA
  )
  
  expect_identical(length(bl), 20L)
  
})

test_that("Parse large covariance matrices", {
  
  test_line <- c("   1.23  4.56   7.89E-02",
                " 1 2 3",
                " 4")
  
  non_square_test <- rep(test_line, 2)
  
  expect_error(
    parse.txtmat(non_square_test),
    "square matrix"
  )
  expect_error(
    ns_mat <- parse.txtmat(non_square_test, square = FALSE),
    NA
  )
  expect_identical(dim(ns_mat), c(2L, 7L))

  
  square_test <- rep(test_line, 7)
  
  expect_error(
    s_mat <- parse.txtmat(square_test),
    NA
  )
  expect_identical(dim(s_mat), c(7L, 7L))

  
})
