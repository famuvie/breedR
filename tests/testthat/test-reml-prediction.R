### For testing prediction, we perform a cross-validation excercise ###
old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

data(m1)

# Sample individuals and generate a secondary response
# with removed observations
Nobs <- nrow(as.data.frame(m1))
sel.idx <- sample(as.data.frame(m1)$self, Nobs/6)

m1$Data$y <- m1$Data$phe_X
m1$Data$y[sel.idx] <- NA
dat <- as.data.frame(m1)

# Fit with no missings
res.full <- remlf90(fixed   = phe_X ~ sex, 
                    genetic = list(model = 'add_animal', 
                                   pedigree = get_pedigree(m1),
                                   id = 'self'), 
                    spatial = list(model = 'AR', 
                                   coord = coordinates(m1),
                                   rho = c(.9, .9)), 
                    data = dat,
                    method = 'ai')

# Fit with missings
res.pred <- remlf90(fixed   = y ~ sex, 
                    genetic = list(model = 'add_animal', 
                                   pedigree = get_pedigree(m1),
                                   id = 'self'), 
                    spatial = list(model = 'AR', 
                                   coord = coordinates(m1),
                                   rho = c(.9, .9)), 
                    data = dat,
                    method = 'ai')


#### Context: Prediction and cross-validation ####
context("Prediction")

# For a cross-validation, the results will not be exactly equal
# but similar up to a given tolerance
tol = 2e-01

# Variance components
test_that("Estimated variance components are similar", {
  expect_equal(res.full$var, res.pred$var, tolerance = tol)
})

# Prediction of the spatial effect
# qplot(ranef(res.full)$spatial,
#       ranef(res.pred)$spatial) +
#   geom_abline(int=0, sl=1)
test_that("Predicted spatial effects are similar everywhere", {
  expect_equal(ranef(res.full)$spatial,
               ranef(res.pred)$spatial,
               tolerance = tol, check.attributes = FALSE)
})



# Prediction of the genetic effect
# plotdat <- data.frame(fullBV = ranef(res.full)$genetic,
#                       predBV = ranef(res.pred)$genetic,
#                       miss   = factor(0, levels = 0:1))
# plotdat$miss[sel.idx] <- 1L
# qplot(fullBV, predBV, data = plotdat, color = miss) + 
#   geom_abline(int = 0, sl = 1, col ='darkgray')
test_that("Predicted Breeding Values are similar in observed individuals", {
  expect_equal(ranef(res.full)$genetic[-sel.idx],
            ranef(res.pred)$genetic[-sel.idx],
            tolerance = tol)
})

test_that("Predicted Breeding Values of missings are still similar 
          up to one additional order of magnitude", {
  expect_equal(ranef(res.full)$genetic[sel.idx],
               ranef(res.pred)$genetic[sel.idx],
               tolerance = 10*tol)
})

# Estimation of fixed effects: can change up to one additional order of magnitude
test_that("Fixed effects are similar", {
  expect_equal(fixef(res.full), fixef(res.pred), tolerance = 10*tol)
})

