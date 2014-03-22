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

# considered as random effects, the factors f3 and o4 have variances
sigma2_f3 <- var(f3.val)
sigma2_o4 <- var(o4.val)

#   # I recover the right coefficients with a standard lm
#   lmres <- lm(y ~ 0 + x + f3 + o4, data = dat)
#   summary(lmres)



#### Context: Linear Models ####
context("Linear Models") 

lm_models <- list(y ~ f3,
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
  # TODO: I can compare variance componens as well.
  
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
  res_em.f3_char <- remlf90(y ~ f3, 
                         data = within(dat, f3 <- as.character(f3)), 
                         method = 'em')
  res_ai.f3_char <- remlf90(y ~ f3, 
                         data = within(dat, f3 <- as.character(f3)), 
                         method = 'ai')
  
  expect_that(is.factor(res_em.f3_char$mf$f3), is_true())
  expect_that(is.factor(res_ai.f3_char$mf$f3), is_true())
})

test_that("remlf90() estimates matches lm()'s", {
  lapply(lm_models, run_expectations, method = 'em')
})

test_that("airemlf90() estimates matches lm()'s", {
  lapply(lm_models, run_expectations, method = 'ai')
})

# Standard Errors
# Notice that while the MLE estimates of beta are equal for continuous
# variables, the standard errors given by aireml and lm are quite different:
# 9.54 and 0.1048 respectively


#### Context: Linear Mixed Models ####
context("Linear Mixed Models") 
require(lme4)

# Re-use lm_models terms as fixed and the available
# factors as random
lm2lmm <- function(x) {
  fixed = x
  random = NULL
  isF <- sapply(dat, is.factor)
  freeF <- names(isF)[which(isF)]
  usedF <- match(attr(terms(x), 'term.labels'), names(isF)[which(isF)])
  usedF <- usedF[!is.na(usedF)]
  if(length(usedF)) freeF <- freeF[-usedF]
  if(length(freeF))
    random = as.formula(paste('~', paste(freeF, collapse = '+')),
                        env = parent.frame(2))
  return(list(fixed = fixed, random = random))
}
lmm_models <- lapply(lm_models, lm2lmm)

# Don't run again models with NULL random terms
lmm_models <- lmm_models[!sapply(lmm_models, function(x) is.null(x$random))]

# Run REML and lm and save estimates and MLEs
run_lmm <- function(m, data = dat, method) {
  res.reml <- remlf90(fixed = m$fixed, 
                      random = m$random, 
                      data = data, 
                      method = method)
  fml.lme4 <- lme4_fml(m$fixed, m$random)
  res.lmm   <- lmer(fml.lme4, data = data)
  return(list(res.reml, res.lmm))
}

# Compare progsf90 and lme4 results
run_lmmexpectations <- function(m, data = dat, method, tol = 1e-03) {
  res <- run_lmm(m, data, method)
  
  # progsf90 estimates should match exactly the
  # MLE estimates from lmer() with the intercept removed.
  
  # TODO: compare standard errors of effects and variance components.
  
  #     # equal standard errors (with more tolerance)
  #     pf90.se <- coef(res[[1]])[[1]]$'s.e.'
  #     lm.se   <- coef(summary(res[[2]]))[, 'Std. Error']
  #     expect_that(pf90.se,
  #                 equals(lm.se, check.attributes = FALSE, tolerance = 1e-05))
  #     
  
  # equal fixed and random effects estimates
  # note that the order of the components might differ
  pf90.beta <- do.call(rbind, coef(res[[1]]))$value
  ranef_order <- names(ranef(res[[1]]))
  lmm.beta   <- c(fixef(res[[2]]), 
                  do.call(rbind, 
                          ranef(res[[2]])[ranef_order])[['(Intercept)']])
  expect_that(pf90.beta,
              equals(lmm.beta, check.attributes = FALSE, tolerance = tol))
  
  # equal variance components estimations
  # note that the order of the components might differ
  pf90.var <- res[[1]]$var[, 1]
  # exclude the residual term, which is extracted separately from lmer
  varnames_order <- head(names(pf90.var), n = -1L)
  lmm.var  <- VarCorr(res[[2]])
  lmm.var  <- c(sapply(lmm.var, function(x) as.numeric(x))[varnames_order],
                Residual = attr(lmm.var, 'sc')^2)
  expect_that(as.numeric(pf90.var),
              equals(lmm.var, check.attributes = FALSE, tolerance = tol))
  
  # equal fitted values
  pf90.fitted <- fitted(res[[1]])
  lmm.fitted   <- fitted(res[[2]])
  expect_that(pf90.fitted, 
              equals(lmm.fitted, check.attributes = FALSE, tolerance = tol))
  qplot(pf90.fitted, lmm.fitted) + geom_abline(int = 0, sl = 1)
}

# Run character conversion test for each method
test_that("Character variables in random effects are treated as factors", {
  res_em.f3_char <- remlf90(fixed = y ~ x, random = ~ f3,
                         data = within(dat, f3 <- as.character(f3)), 
                         method = 'em')
  res_ai.f3_char <- remlf90(fixed = y ~ x, random = ~ f3,
                         data = within(dat, f3 <- as.character(f3)), 
                         method = 'em')
  
  expect_that(is.factor(res_em.f3_char$mf$f3), is_true())
  expect_that(is.factor(res_ai.f3_char$mf$f3), is_true())
})

# Run all test for each of the models and each method
for( m in lmm_models ) {
  m_str <- paste(lapply(m, deparse), collapse = ' | ')

  test_that(paste("remlf90() estimates matches lmer()'s for model", m_str), {
    run_lmmexpectations(m, data = dat, method = 'em')
  })
  
  test_that(paste("airemlf90() estimates matches lmer()'s for model", m_str), {
    run_lmmexpectations(m, data = dat, method = 'ai')
  })
  
}

