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

test_that("The isotropic empirical variograms are well computed", {
  expect_equal(as.numeric(res$variogram), res$truev, tolerance = .2)
})



test_that("The anisotropic empirical variograms are well computed", {
  
  ## Toy test data:
  ## Non-rectangular; not equally-spaced; non-integer coords
  #
  #       522.0| ·  ·  ·  7
  #       518.5| ·  ·  5  6
  #       515.0| 1  2  3  4
  #             ___________
  #             18 20 22 24
        
  coord <- data.frame(x = c(18 + 2*(0:3), c(22, 24, 24)),
                     y = c(rep(515, 4), rep(518.5, 2), 522))
  z <- 1:7
  R <- ceiling(max(dist(coord)))  # compute all combinations
  # ggplot(tdat, aes(x, y, label = z)) + geom_text()
  tvar <- variogram(coord = coord, z = z, R = R)
  
  ## Manually computed results
  ord <- with(tvar$anisotropic, order(x, y))
  true.var <- data.frame(
    N = c(0, 0, 0, 0, 1, 0, 3, 1, 4, 3, 1, 2, 2, 1, 1, 1, 1), 
    v = c(NaN, NaN, NaN, NaN, 1, NaN, 3, 9, 1, 22/3, 16, 4, 16, 25, 9, 25, 36)/2)
  
  expect_equal(tvar$anisotropic[ord, 'n'], true.var$N)
  expect_equal(tvar$anisotropic[ord, 'z'], true.var$v)
})


