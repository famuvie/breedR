old.op <- options(warn = -1,  # suppressWarnings
                  show.error.messages = FALSE)  # silent try
on.exit(options(old.op))

context("Blocks infrastructure")
########################

test_that("breedr_blocks() constructor gives a list with 5 elements of correct sizes", {
  x.loc <- 1:100
  y.loc <- seq(1000, by = 5, length = 51)
  n.blocks <- 10
  
  full.coord <- expand.grid(x.loc, y.loc)
  coord <- full.coord[sample(nrow(full.coord), nrow(full.coord)/2), ]
  block <- factor(sample(n.blocks, nrow(coord), replace = TRUE))

  result <- breedr_blocks(coord, id = block)
  inc.mat <- model.matrix(result)
  cov.mat <- get_structure(result)

  expect_is(result, c('blocks', 'spatial, random', 'breedr_effect'))
  expect_equal(length(result), 5)
  expect_equal(names(result$param), 'n.blocks')
  expect_equal(n.blocks, result$param$n.blocks)
  expect_equal(nrow(inc.mat), nrow(coord))
  expect_equal(ncol(inc.mat), n.blocks)
  
  # The matrix U should be in sparse format: row col value
  expect_that(cov.mat, is_a('Matrix'))
})


context("Extraction of results from spatial blocks model")
########################

data(globulus)
dat <- globulus

fixed.fml <- phe_X ~ gg

n.obs     <- nrow(dat)
n.fixed   <- length(attr(terms(fixed.fml), 'term.labels'))
nlevels.fixed <- nlevels(dat$gg)
n.blocks   <- nlevels(dat$bl)

res <- try(
  suppressMessages(
    remlf90(
      fixed = fixed.fml, 
      spatial = list(model = 'blocks', 
                     coord = globulus[, c('x', 'y')],
                     id = dat$bl), 
      data = dat)
  )
)

# # Manual verification of block estimates:
# library(dplyr)
# fixef(res)
# globulus %>% group_by(gg) %>% summarise(group_mean = mean(phe_X))

# Debug
# tb <- breedr_blocks(globulus[, c('x', 'y')], dat$bl)
# 
# effpf90 <- renderpf90.breedr_modelframe(res$effects, 1)
# pf90 <- progsf90(res$mf, res$effects, opt = '', res.var.ini = 10)

test_that("The blocks model runs with EM-REML without errors", {
  expect_that(!inherits(res, "try-error"), is_true())
})

test_that("coef() gets a named vector of coefficients", {
  expect_is(coef(res), 'numeric')
  expect_equal(length(coef(res)), nlevels.fixed + n.blocks)
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
  expect_named(x, names(res$effects))
  expect_equal(dim(x$gg), c(n.obs, nlevels.fixed))
  expect_is(x$spatial, 'sparseMatrix')
  expect_equal(dim(x$spatial), c(n.obs, n.blocks))
})


test_that("nobs() gets the number of observations", {
  expect_equal(nobs(res), n.obs)
})

test_that("plot(, type = *) returns ggplot objects", {
  expect_is(plot(res, type = 'phenotype'), 'ggplot')
  expect_is(plot(res, type = 'fitted'), 'ggplot')
  expect_is(plot(res, type = 'spatial'), 'ggplot')
  expect_is(plot(res, type = 'fullspatial'), 'ggplot')
  expect_is(plot(res, type = 'residuals'), 'ggplot')
})

test_that("print() shows some basic information", {
  ## Not very informative currently...
  expect_output(print(res), 'Data')
})

test_that("ranef() gets a ranef.breedR object with random effect BLUPs and their s.e.", {
  x <- ranef(res)
  expect_is(x, 'ranef.breedR')
  expect_equal(length(x), 1)
  expect_named(x, c('spatial'))

  expect_is(x$spatial, 'numeric')
  expect_equal(length(x$spatial), n.blocks)
  expect_false(is.null(xse <- attr(x$spatial, 'se')))
  
  expect_is(xse, 'numeric')
  expect_equal(length(xse), n.blocks)
})

test_that("residuals() gets a vector of length N", {
  rsd <- residuals(res)
  expect_is(rsd, 'numeric')
  expect_equal(length(rsd), n.obs)
})

test_that("summary() shows summary information", {
  expect_output(summary(res), 'Variance components')
  expect_output(summary(res), 'blocks')
})

test_that("vcov() gets the covariance matrix of the spatial component of the observations", {
  x <- vcov(res)
  expect_is(x, 'Matrix')
  expect_equal(dim(x), rep(n.obs, 2))
})

