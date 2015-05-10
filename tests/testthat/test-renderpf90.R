old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

context("renderpf90")

test.inc8x4 <- as(rep(1:4, 2), 'indMatrix')   # An 8x4 incidence matrix
test.cov4x4 <- with(L <- Matrix::tril(matrix(sample(16),4)),
                    Matrix::t(L)%*%L)         # a SPD 4x4 covariance matrix

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


## TODO: Test with objects of types: 
##   - competition
##   - genetic
##   - permanent_environmental_competition
##   - additive_genetic
##   - additive_genetic_animal
##   - additive_genetic_competition
