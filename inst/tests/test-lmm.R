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
               y ~ o4)
# Run REML and lm and save estimates and MLEs
run_model <- function(m, data = dat, method) {
  res.reml <- remlf90(m, data = data, method = method)
  res.lm   <- lm(update(m, ~ . -1), data = data)
  return(list(fixef(res.reml)[[1]]$value,
              as.numeric(coefficients(res.lm))))
}

# Compare progsf90 estimates MLEs from lm
compare_lm <- function(m, data = dat, method) {
  res <- run_model(m, data, method)
  expect_that(res[[1]], equals(res[[2]]))
}

test_that("Character variables are treated as factors", {
  res.f3_char <- remlf90(y ~ f3, 
                         data = within(dat, f3 <- as.character(f3)), 
                         method = 'em')
  
  expect_that(is.factor(res.f3_char$mf$f3), is_true())
})

test_that("remlf90 estimates mathes lm's", {
  lapply(models, compare_lm, method = 'em')
})


# Standard Errors
# Notice that while the MLE estimates of beta are equal for continuous
# variables, the standard errors given by reml and lm are quite different:
# 9.54 and 0.1048 respectively

test_that("airemlf90estimates mathes lm's", {
  lapply(models, compare_lm, method = 'ai')
})