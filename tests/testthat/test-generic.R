old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

### Test the generic model ###

context("Generic model")

test.inc8x4 <- as(rep(1:4, 2), 'indMatrix')   # An 8x4 incidence matrix
test.cov4x4 <- with(L <- Matrix::tril(matrix(sample(16),4)),
                    Matrix::t(L)%*%L)         # a SPD 4x4 covariance matrix

test_that('generic() takes an incidence matrix in all matrix-like formats', {
  
  ## logical matrix (converted to Matrix::lgCMatrix)
  logical.matrix <- as.matrix(test.inc8x4)
  gm1 <- generic(incidence  = logical.matrix, covariance = test.cov4x4)
  expect_identical(gm1$incidence.matrix, as(logical.matrix, 'Matrix'))
  
  ## numeric dataframe (convert to matrix and then to Matrix::dgCMatrix)
  numeric.dataframe <- with(list(logical.matrix), {
    mode(logical.matrix) <- 'double'
    data.frame(logical.matrix)})
  gm2 <- generic(incidence  = numeric.dataframe, covariance = test.cov4x4)
  expect_identical(gm2$incidence.matrix,
                   as(as.matrix(numeric.dataframe), 'Matrix'))
  
  ## index Matrix (not converted at all)
  gm3 <- generic(incidence  = test.inc8x4, covariance = test.cov4x4)
  expect_identical(gm3$incidence.matrix, test.inc8x4)
  
})


test_that('generic() takes either covariance or precision matrices', {
  
  gm21 <- generic(incidence  = test.inc8x4, covariance = test.cov4x4)
  expect_match(gm21$structure.type, 'covariance')
  
  gm22 <- generic(incidence  = test.inc8x4, precision = test.cov4x4)
  expect_match(gm22$structure.type, 'precision')
  
})

## TODO:
## Test more throughly:
#   - all.equal.remlf90()
#   - breedr_effect()
#   - effect_group()
#   - model.matrix.effect_group()
#   - model.matrix.breedr_effect()
#   - random()
#   - vcov.random()




context("Extraction of results from generic model")
########################

dat <- globulus

fixed.fml <- phe_X ~ gg + x
n.obs     <- nrow(dat)
n.fixed   <- length(attr(terms(fixed.fml), 'term.labels'))
nlevels.fixed <- nlevels(dat$gg) + 1
nlevels.random <- nlevels(dat$bl)

inc.mat <- model.matrix(~ 0 + bl, globulus)
cov.mat <- diag(nlevels(globulus$bl))

res <- try(remlf90(fixed   = fixed.fml,
                   generic = list(bl = list(inc.mat,
                                            cov.mat)),
                   data    = dat),
           silent = TRUE)



test_that("The generic model runs with AI-REML without errors", {
  expect_that(!inherits(res, "try-error"), is_true())
})

test_that("coef() gets a named vector of coefficients", {
  expect_is(coef(res), 'numeric')
  expect_equal(length(coef(res)), nlevels.fixed + nlevels.random)
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

test_that("get_pedigree() returns NULL", {
  expect_null(get_pedigree(res))
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
  expect_named(x$random, c('bl'))
  expect_is(x$random$bl, 'sparseMatrix')
  expect_equal(dim(x$random$bl), c(n.obs, nlevels.random))
})

test_that("nobs() gets the number of observations", {
  expect_equal(nobs(res), n.obs)
})

test_that("plot(, type = *) returns ggplot objects after providing coords", {
  ## An error mesage is expected as the spatial structure is missing
  expect_error(suppressMessages(plot(res, type = 'phenotype')),
               'Missing spatial structure')
  
  ## We can still plot phenotype, fitted and residuals if provide coords
  coordinates(res) <- dat[, c('x', 'y')]
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
  expect_named(x, c('bl'))
  
  expect_is(x$bl, 'numeric')
  expect_equal(length(x$bl), nlevels.random)
  expect_false(is.null(xse <- attr(x$bl, 'se')))
  
  expect_is(xse, 'numeric')
  expect_equal(length(xse), nlevels.random)
})

test_that("residuals() gets a vector of length N", {
  rsd <- residuals(res)
  expect_is(rsd, 'numeric')
  expect_equal(length(rsd), n.obs)
})

test_that("summary() shows summary information", {
  expect_output(summary(res), 'Variance components')
})

test_that("vcov() gets the covariance matrix of the bl component of the observations", {
  
  ## Make it available after refactoring
  ## when we can recover the structure and model matrices
  expect_error(x <- vcov(res, effect = 'bl'), 'Error in match.arg')
  #   expect_is(x, 'Matrix')
  #   expect_equal(dim(x), rep(n.obs, 2))
})

