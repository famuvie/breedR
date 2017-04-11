old.op <- options(warn = -1,  # suppressWarnings
                  show.error.messages = FALSE)  # silent try
on.exit(options(old.op))

context("Interface for weights")
########################

test_that("simulated example of WEIGHTS option", {
  
  ## simulated dataset
  set.seed(123)
  n <- 1e4
  res_var <- 4
  dat <- 
    transform(
      data.frame(
        w = runif(n, min = .5, max = 2),
        e = rnorm(n, sd = sqrt(res_var))
      ),
      y = 10 + e/sqrt(w)
    )
  
  ## Now, if R = V(e),
  ## the residual variance for y is R/w
  res <- remlf90(
    y ~ 1,
    data = dat,
    weights = dat$w
  )
  
  Rest <- res$var[1, 1]
  
  # summary(res)
  # plot(dat$e, residuals(res), pch = ".")
  # abline(0, 1)

  ## The sample variance estimator s^2 is such that
  ## (n-1)*s^2/\sigma^2 ~ Chi^2(n-1)
  ## Expect that the estimation lies in the centered 1% confidence interval
  expect_within <- function(value, interval) {
    expect_true(all( c(1, -1) * (value - interval) > 0 ))
  }
  expect_within((n-1)*Rest/res_var, qchisq(c(.005, .995), df = n-1))
})
