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
o4 <- factor(sample(month.abb[1:4], n, replace = TRUE), 
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


test_that("remlf90 fits models with one single covariate", {
  # response against one factor
  res.f3 <- remlf90(y ~ f3, data = dat, method = 'em')
  expect_that(fixef(res.f3)$f3$value,
              is_equivalent_to(as.numeric(by(dat$y, dat$f3, mean))))

  # Character variables are treated as factors
  res.f3_char <- remlf90(y ~ f3, 
                         data = within(dat, f3 <- as.character(f3)), 
                         method = 'em')
  expect_that(res.f3[-1], equals(res.f3_char[-1]))
  
#   # Issue regarding the number of levels of covariates
#   res.x <- remlf90(y ~ x, data = dat)
#   expect_that(res.x,
#               equals(test.length))
})

test_that("airemlf90 gives the same results as remlf90", {
  
})