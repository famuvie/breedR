context("Basic Multitrait models")

## Bivariate model with a random effect of block and a genetic effect.
res_mt <- readRDS(file.path(testdata, "res_mt.rds"))
ntraits <- ncol(as.matrix(model.response(res_mt$mf)))

test_that("Extract random effects", {
  
  ## ranef() recovers a list of all existing random effects
  expect_identical(length(ranef(res_mt)), length(res_mt$ranef))
  
  ## each element is a matrix with as many columns as traits
  expect_identical(lapply(ranef(res_mt), ncol), lapply(res_mt$ranef, length))
})



context("Multitrait-competition models")

## Bivariate model with a random effect of block and a genetic effect.
res_mtcp <- readRDS(file.path(testdata, "res_mtcp.rds"))
ntraits <- ncol(as.matrix(model.response(res_mtcp$mf)))

test_that("Extract random effects", {
  
  ## ranef() recovers a list of all existing random effects
  expect_identical(length(ranef(res_mtcp)), length(res_mtcp$ranef))
  
  ## each element is a matrix with as many columns as traits
  expect_identical(lapply(ranef(res_mtcp), ncol), lapply(res_mtcp$ranef, length))
})



