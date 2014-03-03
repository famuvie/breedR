context("Number of knots")
##########################

test_that("determine.n.knots works for atomic vectors", {
  test.length <- 100
  sample.sizes <- seq(from = 7, by = 19, length = test.length)
  expect_that(length(breedR:::determine.n.knots(sample.sizes)),
              equals(test.length))
})

test_that("determine.n.knots fails with few data points", {
  expect_that(breedR:::determine.n.knots(6), throws_error('few data'))
})



context("Splines infraestructure") 
########################

test_that("build.splines.model gives a list with tree elements of correct sizes", {
  x.loc <- 1:100
  y.loc <- seq(1000, by = 5, length = 51)
  coord <- expand.grid(x.loc, y.loc)
  result <- breedR:::build.splines.model(coord)
  n.knots <- ncol(result$B)
  
  expect_that(result, is_a('list'))
  expect_that(length(result), equals(5))
  expect_that(nrow(result$B), equals(nrow(coord)))
  expect_that(prod(result$inner.knots + 2), equals(n.knots))
  # The matrix U should be in sparse format: row col value
  expect_that(ncol(result$U), equals(3))
  expect_that(max(result$U[, 1]), equals(n.knots))
  expect_that(max(result$U[, 2]), equals(n.knots))
  
  # The last element is a list of coordinates
  # and the corresponding incidence matrix
  # for plotting purposes
  expect_that(result$plotting, is_a('list'))
  expect_that(dim(result$plotting$B), 
              equals(c(nrow(result$plotting$grid), n.knots)))
})


context("Spatial splines model")
########################

data(m1)
dat <- as.data.frame(m1)

res <- try(
  remlf90(fixed = phe_X ~ sex, 
          spatial = list(model = 'Cappa07', 
                         coord = coordinates(m1),
                         var.ini = 300,
                         n.knots = c(2, 2)), 
          data = dat,
          method = 'em'),
  silent = TRUE
)


test_that("The Cappa07 model runs with EM-REML without errors", {
  expect_that(!inherits(res, "try-error"), is_true())
})