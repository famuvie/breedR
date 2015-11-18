old.op <- options(warn = -1,  # suppressWarnings
                  show.error.messages = FALSE)  # silent try
on.exit(options(old.op))


#### Context: breedr_ar() ####
context("AR infrastructure")

# Some coordinates, with a whole row and some points missing
# with non-integer and non-positive value
test.pos <- list(x = c(-2:3), y = 5:8/2)
test.coord <- as.matrix(expand.grid(test.pos$x[-2], test.pos$y))[-8, ]
  # plot(test.coord, pch = 19)
test.rho <- c(.6, .8)
reslst <- list(breedr_ar(test.coord, test.rho, TRUE),
               breedr_ar(test.coord, test.rho, FALSE))

check_build.ar.model <- function(x) {
  inc.mat <- model.matrix(x)
  cov.mat <- get_structure(x)
  eff.size <- ifelse(attr(x, 'grid')$autofill,
                     prod(sapply(test.pos, length)),
                     nrow(test.coord)+1)
  
  test_that("breedr_ar() returns a list with the right elements", {
    expect_is(x, c('ar', 'spatial', 'random', 'breedr_effect'))
    expect_equal(length(x), 5)
    expect_equal(names(x$param), 'rho')
    expect_equal(x$param$rho, test.rho)
    # Incidence matrix
    expect_is(inc.mat, 'sparseMatrix') # a permutation Matrix
    expect_equal(nrow(inc.mat), nrow(test.coord))
    # Covariance matrix
    expect_is(cov.mat, 'sparseMatrix') 
    expect_equal(ncol(cov.mat), eff.size)
    expect_equal(ncol(inc.mat), nrow(cov.mat))
  })
}

for (x in reslst) check_build.ar.model(x)



#### Context: AR models with different arrangements of trees ####
context("AR models with diffferent arrangements of trees")


#### Build small testbeds ####
build.testbed <- function(corner = c(0, 0), size, treesep = c(1, 1), beta){
  n = size[1] * size[2]
  # A planar spatial effect
  s.mat = matrix(NA, nrow = size[1], ncol = size[2])
  j = 1:size[2]
  for(i in 1:size[1])
    s.mat[i,j] = 0.1*(i+2*j)
  
  ## a covariate
  set.seed(2)
  z.mat = matrix(runif(n), size[1], size[2])
  ## noise
  set.seed(2)
  noise.mat = matrix(rnorm(n, sd = 0.3), size[1], size[2])
  ## make simulated data
  y.mat = beta * z.mat + s.mat + noise.mat
  ## build final dataset  
  dat <- data.frame(i = rep(seq(corner[1], by = treesep[1], length = size[1]),
                            times = size[2]),
                    j = rep(seq(corner[2], by = treesep[2], length = size[2]),
                            each = size[1]),
                    z = as.vector(z.mat),
                    y = as.vector(y.mat),
                    true.s = as.vector(s.mat))
  return(dat)
}

beta = 0.5
datlist <- list(# small square regular grid
  small.sq.reg = build.testbed(corner = c(0, 0),
                               size = c(5, 5),
                               treesep = c(1, 1),
                               beta = beta),
  # small rectangular grid with different spacings and coordinates
  small.rect.irr = build.testbed(corner = c(134, 77),
                                 size = c(5, 7),
                                 treesep = c(3, 4),
                                 beta = beta))

# triangular configuration
datlist <- c(datlist,
             triang = list(datlist[[1]][which(as.vector(Matrix::tril(matrix(TRUE, 5, 5)))),]))




# Fit models both with EM and AI-REML
run.model <- function(dat, method) {
  res = try(
    suppressMessages(
      remlf90(
        fixed = y ~ 1 + z, 
        spatial = list(model = 'AR',
                       coord = dat[, 1:2],
                       rho = c(.9, .9)),
        data = dat,
        method = method)
    )
  )
  return(list(dat = dat,
              method = method,
              res = res))
}

reslist <- c(lapply(datlist, run.model, method = 'em'),
             lapply(datlist, run.model, method = 'ai'))
             
# Check results
# summary(reslist[[1]])
# res <- reslist[[1]]
# dat <- datlist[[1]]
# require(plyr)
check.result <- function(m, datlabel, debug.plot = FALSE) {
  test_that(paste("AR model runs OK with dataset", datlabel, "and method", m$method), {
    expect_true(!inherits(m$res, 'try-error'))
  })
  
  if( !inherits(m$res, 'try-error') ){
    fit.s <- fixef(m$res)$Intercept$value +
      model.matrix(m$res)$spatial %*% ranef(m$res)$spatial
    if(debug.plot) {
      print(qplot(as.vector(m$dat$true.s), fit.s) +
              geom_abline(int = 0, sl = 1))
    }
    # Mean Square Error for the spatial effect
    mse <- mean((as.vector(m$dat$true.s) - fit.s)^2)
    test_that(paste("MSE of the spatial effect estimation is reasonable for dataset",
                    datlabel, "and method", m$method), {
      expect_that(mse, is_less_than(1))
    })
    
    # Estimate of the linear coefficient
    beta.e <- beta - fixef(m$res)$z$value
    test_that(paste("The linear coefficient is estimated within 3 se for dataset",
              datlabel, "and method", m$method), {
      expect_that(abs(beta.e), is_less_than(3*fixef(m$res)$z$s.e.))
    })
  }
}

for(i in 1:length(reslist)) 
  check.result(reslist[[i]], names(reslist)[i], debug.plot = FALSE)



#### Context: selection of autoregressive parameters ####
context("Selection of autoregressive parameters")

res.unset <- try(
  suppressMessages(
    remlf90(
      fixed = y ~ z, 
      spatial = list(model = 'AR',
                     coordinat = datlist[[1]][, 1:2]),
      data = datlist[[1]])
  )
)

test_that("if rho unset, remlf90 tries a grid of combinations", {
  # remlf90() returns an evaluation grid
  expect_that(exists('rho', as.environment(res.unset)), is_true())
  # the evaluation grid returns the loglikelihood for each default combination
  expect_that(all(complete.cases(res.unset$rho$loglik)), is_true())
})


gridlist <- list(expand.grid(seq(80, 90, 5), c(87, 93))/100,
                 expand.grid(seq(80, 90, 5), NA)/100,
                 expand.grid(NA, c(87, 93))/100)
reslist.spec <- lapply(gridlist, function(g)
  try(
    suppressMessages(
      remlf90(
        fixed = y ~ z, 
        spatial = list(model = 'AR',
                       coord = datlist[[1]][, 1:2],
                       rho = g),
        data = datlist[[1]])
    )
  )
)

test_that("the user can specify a full or partial grid of combinations", {
  
  for(i in 1:length(reslist.spec)) {
    res <- reslist.spec[[i]]
    grid <- gridlist[[i]]
    
    # remlf90() returns an evaluation grid
    expect_that(exists('rho', as.environment(res)), is_true())
    
    # The evaluation grid conforms to the user specification
    get_levels <- function(levels) {
      if(all(is.na(levels))) return(breedR.getOption('ar.eval'))
      else return(levels)
    }
    eval.grid <- expand.grid(lapply(lapply(grid, unique), get_levels),
                             KEEP.OUT.ATTRS = FALSE)
    names(eval.grid) <- names(res$rho)[1:2]
    expect_identical(eval.grid,
                     res$rho[, 1:2])

    # the evaluation grid returns the loglikelihood for each combination specified
    expect_that(all(complete.cases(res$rho$loglik)),
                is_true())
  }
})


# # Debug
# image(s.mat)
# image(matrix(res.bR$spatial$fit$z, nrow, ncol))
# qplot(as.vector(s.mat), res.bR$spatial$fit$z) + geom_abline(int = 0, sl = 1, col = 'darkgray')
# summary(res.bR)




context("Extraction of results from spatial AR model")
########################

data(m1)
dat <- as.data.frame(m1)

## Remove some observations to provoke
## misalignment beetween the observations and the spatial random effects
dat <- dat[-sample(1:nrow(dat), 50), ]

fixed.fml <- phe_X ~ sex

n.obs     <- nrow(dat)
n.fixed   <- length(attr(terms(fixed.fml), 'term.labels'))
nlevels.fixed <- nlevels(dat$sex)
rho   <- c(.9, .9)
coord = dat[, 1:2]

## Number of levels of the AR effect
n.AR <- prod(sapply(loc_grid(coord, autofill = TRUE), length))

# Use a different number of knots for rows and columns
res <- try(
  suppressMessages(
    remlf90(fixed = fixed.fml, 
            spatial = list(model = 'AR', 
                           coord = coord,
                           rho = rho), 
            data = dat,
            method = 'ai')
  )
)


test_that("The AR model runs with EM-REML without errors", {
  expect_that(!inherits(res, "try-error"), is_true())
})

test_that("coef() gets a named vector of coefficients", {
  expect_is(coef(res), 'numeric')
  expect_equal(length(coef(res)), nlevels.fixed + n.AR)
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
  expect_equal(dim(x$sex), c(n.obs, nlevels.fixed))
  expect_is(x$spatial, 'sparseMatrix')
  expect_equal(dim(x$spatial), c(n.obs, n.AR))
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
  expect_equal(length(x$spatial), n.AR)
  expect_false(is.null(xse <- attr(x$spatial, 'se')))
  
  expect_is(xse, 'numeric')
  expect_equal(length(xse), n.AR)
})

test_that("residuals() gets a vector of length N", {
  rsd <- residuals(res)
  expect_is(rsd, 'numeric')
  expect_equal(length(rsd), n.obs)
})

test_that("summary() shows summary information", {
  expect_output(summary(res), 'Variance components')
  expect_output(summary(res), 'rho:')
})

test_that("vcov() gets the covariance matrix of the spatial component of the observations", {
  x <- vcov(res)
  expect_is(x, 'Matrix')
  expect_equal(dim(x), rep(n.obs, 2))
})

