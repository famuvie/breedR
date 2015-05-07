old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

context("renderpf90")

##Test with an object of type generic
gtest <- generic(incidence  = test.inc8x4, covariance = test.cov4x4)
test_that("The function returns the same number",{
  expect_equal(renderpf90(gtest)$levels,nrow(gtest$structure.matrix))
})

##Test with an object of type splines
x.loc <- 1:100
y.loc <- seq(1000, by = 5, length = 51)
coord <- expand.grid(x.loc, y.loc)
result <- splines(coord)
test_that("The function behaves as expected",{
  expect_equal(tail(renderpf90(result)$levels, 1),nrow(result$structure.matrix))
  expect_true(all(head(renderpf90(result)$levels, -1) ==0))
})

