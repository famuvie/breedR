old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

#### Context: Simulation infrastructure ####
context("Simulation infrastructure") 


dat <- try(breedR.sample.phenotype(fixed   = c(mu = 10, x = 2),
                                   random = list(u = list(nlevels = 3,
                                                          sigma2  = 1)),
                                   genetic = list(model    = 'competition',
                                                  Nparents = c(10, 10),
                                                  sigma2_a = matrix(c(2, -1, -1, 2), 2, 2),
                                                  competition_decay = 1,
                                                  check.factorial = FALSE,
                                                  pec = 0.5),
                                   spatial = list(model     = 'AR',
                                                  grid.size = c(10, 10),
                                                  rho   = c(.9, .5),
                                                  sigma2_s  = 1),
                                   residual.variance = 1))

test_that('breedR.sample.phenotype() runs without error', {
  expect_false(inherits(dat, 'try-error'))
})