old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

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

test_that("splines() constructor gives a list with six elements of correct sizes", {
  x.loc <- 1:100
  y.loc <- seq(1000, by = 5, length = 51)
  coord <- expand.grid(x.loc, y.loc)
  result <- splines(coord)
  inc.mat <- model.matrix(result)
  cov.mat <- get_structure(result)
  n.knots <- sapply(result$knots, length)
  n.splines <- ncol(cov.mat)
  
  expect_that(result, is_a('splines'))
  expect_that(length(result), equals(6))
  expect_equal(n.knots,
               sapply(sapply(list(x.loc, y.loc), length), determine.n.knots)+6,
               check.attributes = FALSE)
  expect_that(nrow(inc.mat), equals(nrow(coord)))
  expect_that(prod(n.knots-4), equals(n.splines))

  # The matrix U should be in sparse format: row col value
  expect_that(cov.mat, is_a('Matrix'))
  
})


context("Spatial splines model")
########################

data(m1)
dat <- as.data.frame(m1)

# Use a different number of knots for rows and columns
res <- try(remlf90(fixed = phe_X ~ sex, 
                   spatial = list(model = 'splines', 
                                  coord = coordinates(m1),
                                  n.knots = c(2, 3)), 
                   data = dat,
                   method = 'em'),
           silent = TRUE)


test_that("The splines model runs with EM-REML without errors", {
  expect_that(!inherits(res, "try-error"), is_true())
})