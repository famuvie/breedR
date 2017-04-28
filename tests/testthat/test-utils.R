### Test the auxiliar functions in utils.R ###


context("Auxiliar functions")


test_that("lmat2df() works as expected", {
  
  rnms <- letters[1:3]
  cnms <- c("x", "y")
  label <- "e"
  
  tm <- matrix(1:9, 3, 3, dimnames = rep(list(rnms), 2))
  
  tlm <- structure(rep(list(tm), 2), names = cnms)
  
  tdf <- lmat2df(tlm, label)

  ## labels
  exp_labs <- c("e.a", "e.a_e.b", "e.a_e.c", "e.b", "e.b_e.c", "e.c")
  expect_identical(rownames(tdf), exp_labs)
  
  ## columns
  expect_identical(colnames(tdf), cnms)
  
  ## values
  ## only lower-triangular values
  exp_values <- tm[lower.tri(tm, diag = TRUE)]
  expect_identical(tdf$x, exp_values)
  expect_identical(tdf$y, exp_values)
})
