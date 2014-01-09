#### Simulated dataset ####
# dataset size
n <- 1e3
# a covariate
x <- runif(n)
# a factor with 3 levels and associated values
f3 <- factor(sample(letters[1:3], n, replace = TRUE))
f3.val <- -1:1
# an ordered factor with 4 levels and associated values
# following a linear trend
o4 <- ordered(sample(month.abb[1:4], n, replace = TRUE), 
             levels = month.abb[1:4])
o4.val <- 1:4
# residual variance
sigma2 <- .5
dat <- data.frame(
  y = x + f3.val[f3] + o4.val[o4] + rnorm(n, sd = sqrt(sigma2)),
  x, f3, o4)
#   # I recover the right coefficients with a standard lm
#   lmres <- lm(y ~ 0 + x + f3 + o4, data = dat)
#   summary(lmres)


#### Contexts ####
context("Linear Models") 

models <- list(y ~ f3,
               y ~ x,
               y ~ o4,
               y ~ x + f3,
               y ~ x + o4,
               y ~ f3 + o4,
               y ~ x + f3 + o4)
# Run REML and lm and save estimates and MLEs
run_model <- function(m, data = dat, method) {
  res.reml <- remlf90(m, data = data, method = method)
  res.lm   <- lm(update(m, ~ . -1), data = data)
  return(list(res.reml, res.lm))
}

# Compare progsf90 and lm results
run_expectations <- function(m, data = dat, method) {
  res <- run_model(m, data, method)
  
  # Expectations will depend on the number of covariates included in the model
  # For one single covariate, progsf90 estimates should match exactly the
  # MLE estimates from lm() with the intercept removed.
  # For more than one covariate, the model paremeterizations are different
  # and I can't compare the coefficients directly. Rather, I compare the
  # fitted values.

  n.cov <- length(attr(terms(m), 'term.labels'))
  
  if( n.cov == 1) {
    # equal fixed effects estimates
    pf90.beta <- coef(res[[1]])[[1]]$value
    lm.beta   <- coef(res[[2]])
    expect_that(pf90.beta,
                equals(lm.beta, check.attributes = FALSE))
    
    # equal standard errors (with more tolerance)
    pf90.se <- coef(res[[1]])[[1]]$'s.e.'
    lm.se   <- coef(summary(res[[2]]))[, 'Std. Error']
    expect_that(pf90.se,
                equals(lm.se, check.attributes = FALSE, tolerance = 1e-05))
    
    # (almost) equal residual variance estimates
    pf90.sigmahat <- sqrt(as.numeric(res[[1]]$var[1, 1]))
    lm.sigmahat   <- summary(res[[2]])$sigma
    expect_that(pf90.sigmahat, equals(lm.sigmahat, tolerance = 1e-03))
  } else {
    
    # equal fitted values
    pf90.fitted <- fitted(res[[1]])
    lm.fitted   <- fitted(res[[2]])
    expect_that(pf90.fitted, equals(lm.fitted))
  }
}

test_that("Character variables are treated as factors", {
  res.f3_char <- remlf90(y ~ f3, 
                         data = within(dat, f3 <- as.character(f3)), 
                         method = 'em')
  
  expect_that(is.factor(res.f3_char$mf$f3), is_true())
})

test_that("remlf90() estimates matches lm()'s", {
  lapply(models, run_expectations, method = 'em')
})

test_that("airemlf90() estimates matches lm()'s", {
  lapply(models, run_expectations, method = 'ai')
})

# Standard Errors
# Notice that while the MLE estimates of beta are equal for continuous
# variables, the standard errors given by aireml and lm are quite different:
# 9.54 and 0.1048 respectively
