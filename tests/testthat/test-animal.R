
data(m1)
dat <- as.data.frame(m1)
ped <- get_pedigree(m1)

#### Context: Building additive-genetic models ####
context("Building additive-genetic models")

test_that("Correctly builds structure of additive_genetic_animal component", {
  dat <- data.frame(id = 1:4,
                    sire = c(11, 11, 2, 3),
                    dam  = c(12, NA, 1, 12))
  ## Recoded pedigree with further unobserved descendants
  ## (with will come later in the codification)
  ped <- suppressWarnings(build_pedigree(1:3, data = rbind(dat, c(5, 1, 2))))
  aga <- try(breedR:::additive_genetic_animal(ped, dat$id))
  
  expect_true(!inherits(aga, "try-error"))
  expect_is(aga, 
            c("additive_genetic_animal", 
              "additive_genetic",
              "genetic",
              "random", 
              "breedr_effect"))
  expect_named(aga, 
               c("incidence.matrix", 
                 "structure.matrix",
                 "structure.type", 
                 "pedigree"))
  expect_identical(dim(aga$incidence.matrix),
                   c(nrow(dat), nrow(as.data.frame(ped))))
  expect_identical(dim(aga$structure.matrix),
                   rep(nrow(as.data.frame(ped)), 2))
  expect_identical(aga$structure.type, 'covariance')
  expect_identical(aga$pedigree, ped)
})


#### Context: Animal Models ####
context("Results from Animal Models") 

fixed_models <- list(phe_X ~ sex)

# Run REML and lm and save estimates and MLEs
run_model <- function(m, data = dat, method) {
  res.reml <- try(
    suppressMessages(
      remlf90(fixed = m,
              genetic = list(model = 'add_animal', 
                             pedigree = ped,
                             id = 'self'), 
              data = data,
              method = method)
    )
  )
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
  #  - compare the estimated genetic and residual vaiances with true values
  #  - compare the estimated and true Breeding Values
  #  - compare results to those from package pedigreemm 
  #    (an extension to lme4 to include animal models)
  
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

test_that("model.matrix() gets a named list of incidence matrices", {
  x <- model.matrix(res)
  expect_is(x, 'list')
  expect_named(x, names(res$effects))
  expect_equal(dim(x$sex), c(n.obs, nlevels.fixed))
  expect_is(x$genetic, 'sparseMatrix')
  expect_equal(dim(x$genetic), c(n.obs, n.bvs))
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
  expect_output(print(summary(res)), 'Variance components')
})

test_that("vcov() gets the covariance matrix of the genetic component of the observations", {
  
  x <- try(vcov(res, effect = 'genetic'))
  expect_false(inherits(x, 'try-error'))
  expect_is(x, 'Matrix')
  expect_equal(dim(x), rep(n.obs, 2))
})

