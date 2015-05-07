## Test of the function effect_type 

old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

context("effect_type")

##Test with an object of type generic
gtest <- generic(incidence  = test.inc8x4, covariance = test.cov4x4)
test_that("The function returns the type expected",{
  expect_output(effect_type(gtest),"random")
})

##Test with an object of type splines
x.loc <- 1:100
y.loc <- seq(1000, by = 5, length = 51)
coord <- expand.grid(x.loc, y.loc)
result <- splines(coord)
test_that("The function returns the type expected",{
  expect_output(effect_type(result),"random")
})
