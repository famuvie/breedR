old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

data(m1)
dat <- as.data.frame(m1)
ped <- get_pedigree(m1)

#### Context: Animal Models ####
context("Animal Models") 

fixed_models <- list(phe_X ~ sex)

# Run REML and lm and save estimates and MLEs
run_model <- function(m, data = dat, method) {
  res.reml <- try(
    remlf90(fixed = m,
            genetic = list(model = 'add_animal', 
                           pedigree = ped,
                           id = 'self'), 
            data = data,
            method = method),
    silent = TRUE)
  return(res.reml)
}

# Compare progsf90 and pedigreemm results
run_expectations <- function(m, data = dat, method) {
  res <- run_model(m, data, method)
  
  # It runs without errors 
  test_that("The animal model runs without errors", {
    expect_that(!inherits(res, "try-error"), is_true())
  })
  
  # TODO:
  # other checks, like:
  #  + compare the estimated and true Breeding Values
  #  - compare results to those from package pedigreemm 
  #    (an extension to lme4 to include animal models)

  test_that("The PBVs are around the true values", {
  expect_equal(dat$BV_X, ranef(res)$genetic[-(1:160)], tolerance = 1, check.attributes = FALSE)
})

}



# Run expectations for all models and methods
test_that("remlf90() estimates matches lm()'s", {
  lapply(fixed_models, run_expectations, method = 'em')
})

test_that("airemlf90() estimates matches lm()'s", {
  lapply(fixed_models, run_expectations, method = 'ai')
})



context("Extraction of results from add_animal model")
########################


fixed.fml <- phe_X ~ sex
n.obs     <- nrow(dat)
n.fixed   <- length(attr(terms(fixed.fml), 'term.labels'))
nlevels.fixed <- nlevels(dat$sex)
n.bvs <- nrow(as.data.frame(ped))

res <- run_model(fixed_models[[1]], method = 'ai')


test_that("The add_animal model runs with EM-REML without errors", {
  expect_that(!inherits(res, "try-error"), is_true())
})

test_that("coef() gets a named vector of coefficients", {
  expect_is(coef(res), 'numeric')
  expect_equal(length(coef(res)), nlevels.fixed + n.bvs)
  expect_named(coef(res))
})

test_that("ExtractAIC() gets one number", {
  expect_is(extractAIC(res), 'numeric')
  expect_equal(length(extractAIC(res)), 1)
})

test_that("fitted() gets a vector of length N", {
  expect_is(fitted(res), 'numeric')
  expect_equal(length(fitted(res)), n.obs)
})

test_that("fixef() gets a named list of data.frames with estimated values and s.e.", {
  x <- fixef(res)
  expect_is(x, 'list')
  expect_named(x)
  expect_equal(length(x), n.fixed)
  for (f in x) {
    expect_is(f, 'data.frame')
    expect_named(f, c('value', 's.e.'))
  }
})

test_that("get_pedigree() returns the given pedigree", {
  expect_identical(get_pedigree(res), ped)
})

test_that("logLik() gets an object of class logLik", {
  expect_is(logLik(res), 'logLik')
})

test_that("model.frame() gets an Nx2 data.frame with a 'terms' attribute", {
  x <- model.frame(res)
  expect_is(x, 'data.frame')
  expect_is(terms(x), 'terms')
  expect_equal(dim(x), c(n.obs, n.fixed + 1))
})

test_that("model.matrix() gets a named list of fixed and random incidence matrices", {
  x <- model.matrix(res)
  expect_is(x, 'list')
  expect_named(x, c('fixed', 'random'))
  expect_equal(dim(x$fixed), c(n.obs, nlevels.fixed))
  expect_is(x$random, 'list')
  expect_named(x$random, c('genetic'))
  expect_is(x$random$genetic, 'sparseMatrix')
  expect_equal(dim(x$random$genetic), c(n.obs, n.bvs))
})

test_that("nobs() gets the number of observations", {
  expect_equal(nobs(res), n.obs)
})

test_that("plot(, type = *) returns ggplot objects after providing coords", {
  ## An error mesage is expected as the spatial structure is missing
  expect_error(suppressMessages(plot(res, type = 'phenotype')),
               'Missing spatial structure')

  ## We can still plot phenotype, fitted and residuals if provide coords
  coordinates(res) <- dat[, 1:2]
  expect_is(plot(res, type = 'phenotype'), 'ggplot')
  expect_is(plot(res, type = 'fitted'), 'ggplot')
  expect_is(plot(res, type = 'residuals'), 'ggplot')
  
  ## But still get errors for the absent spatial components
  expect_error(plot(res, type = 'spatial'), 'no spatial effect')
  expect_error(plot(res, type = 'fullspatial'), 'no spatial effect')
})

test_that("print() shows some basic information", {
  ## Not very informative currently...
  expect_output(print(res), 'Data')
})

test_that("ranef() gets a ranef.breedR object with random effect BLUPs and their s.e.", {
  x <- ranef(res)
  expect_is(x, 'ranef.breedR')
  expect_equal(length(x), 1)
  expect_named(x, c('genetic'))
  
  expect_is(x$genetic, 'numeric')
  expect_equal(length(x$genetic), n.bvs)
  expect_false(is.null(xse <- attr(x$genetic, 'se')))
  
  expect_is(xse, 'numeric')
  expect_equal(length(xse), n.bvs)
})

test_that("residuals() gets a vector of length N", {
  rsd <- residuals(res)
  expect_is(rsd, 'numeric')
  expect_equal(length(rsd), n.obs)
})

test_that("summary() shows summary information", {
  expect_output(summary(res), 'Variance components')
})

test_that("vcov() gets the covariance matrix of the genetic component of the observations", {
  
  ## Make it available after refactoring
  ## when we can recover the structure and model matrices
  expect_error(x <- vcov(res, effect = 'genetic'), 'Currently not available')
  #   expect_is(x, 'Matrix')
  #   expect_equal(dim(x), rep(n.obs, 2))
})

