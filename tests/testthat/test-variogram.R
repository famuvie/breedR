old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

### Empirical Variograms ###
context("Empirical Variograms")


### Simulate gaussian observations with a covariance structure produced
### by a given (an)isotropic variogram, and compare the empirical and 
### the theoretical variograms

# Exponential isotropic variogram
# par = list(n = nugget, s = sill, r = range)
truevar.iso <- function(x, par) {
  x[x < 0] <- 0
  gamma <- with(par, (s - n)*(1 - exp(-3*x/r)))
  gamma[x > 0] <- gamma[x > 0] + par$n
  gamma
}

# Covariance function
# 
truecov.iso <- function(x, par) {
  par$s - truevar.iso(x, par)
}

# # check
# par <- list(n = 1, s = 5, r = 10)
# x <- seq(-3, 20, by = .2)
# yv <- truevar.iso(x, par)
# yc <- truecov.iso(x, par)
# plot(x, yv, type = 'l', ylim = c(0, 7))
# lines(x, yc)

# Simulated data
truepar <- list(n = 1, s = 5, r = 20)
coord <- expand.grid(x = 0:49, y = 0:49)
D <- dist(coord)
Slt <- spam::as.spam(truecov.iso(D, par = truepar))
S <- Slt + spam::t(Slt) + spam::diag.spam(truecov.iso(spam::diag(Slt),
                                                      par = truepar))
z <- as.vector(spam::rmvnorm.spam(1, Sigma = S))

# Distances, true and estimated variograms
empvar <- variogram(R = 30, coord = coord, z = z)
res <- transform(empvar$isotropic,
                 truev = truevar.iso(distance, truepar))

# # check
# # Note that the anisotropic variogram is remarkably less robust
# empvar
# with(res, {
#      plot(distance, variogram, type = 'l', ylim = c(0, 7));
#      lines(distance, truev, col ='red')})

test_that("The empirical variograms are well computed", {
  expect_equal(as.numeric(res$variogram), res$truev, tolerance = .1)
})

test_that("variogram needs to be provided with objects inheriting from 'BreedR'",{
  expect_error(variogram(0.5))
  expect_error(variogram(matrix(1:16,4,4)))
  expect_error(variogram(list('test')))
  expect_error(variogram(FALSE))
  })

test_that("if x is not provided in the variogram, then both coord and z must be given'",{
  expect_error(variogram(coord = coord))
  expect_error(variogram(z = z))
})

test_that("variogram runs without error with breedR results or specific coordinates/z",{
  breedtest <- variogram(breedR.result())
  spectest <- variogram(coord = coord, z = z)
  coord2 <- expand.grid(x = 0:29, y = 0 : 29)
  z2 <- rexp(length(coord2$x))
  spectest2 <- variogram(coord = coord2, z = z2)
  expect_that(!inherits(breedtest, "try-error"), is_true())
  expect_that(!inherits(spectest, "try-error"), is_true())
  expect_that(!inherits(spectest2, "try-error"), is_true())
})


