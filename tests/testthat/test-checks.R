old.op <- options(warn = -1,  # suppressWarnings
                  show.error.messages = FALSE)  # silent try
on.exit(options(old.op))

### Tests of functions for checking model components ###

context("Checks for model components")

##Model add_animal##
dat <- data.frame(id = 1:4,
                  sire = c(11, 11, 2, 3),
                  dam  = c(12, NA, 1, 12))
ped <- build_pedigree(1:3, data = dat)
id <- dat$id
var.ini <- 1.5

test_that("The minimal add_animal specification checks without error and completes missing values",{
  
  ## Try alternative correct specification formats:
  add_animal_minimalspec <- list(
    try(check_genetic(model = 'add_animal',
                      pedigree = ped,
                      id = id)),
    try(check_genetic(model = 'add_animal',
                      pedigree = ped,
                      id = 'id',
                      data = dat)),
    try(check_genetic(model = 'add_animal',
                      pedigree = as.data.frame(ped)[-(1:2),],
                      id = id)),
    try(check_genetic(model = 'add',
                      pedigree = ped,
                      id = id))
  )
  
  for (x in add_animal_minimalspec) {
    expect_false(inherits(x, "try-error"))
    
    expect_true(setequal(names(x),
                         c('model', 'pedigree', 'id', 'var.ini', 'autofill')))
    
    ## var.ini should have been added with the default value
    ## and the attribute 'var.ini.default' set to TRUE
    expect_equal(x$var.ini, breedR.getOption('default.initial.variance'))
    expect_true(attr(x, 'var.ini.default'))
  }
})

test_that("check_genetic() returns an error if missing values",{
  expect_error(check_genetic(pedigree = ped, id = id, var.ini = var.ini ))
  expect_error(check_genetic(model = 'add_animal', pedigree = ped, var.ini = var.ini ))
  expect_error(check_genetic(model = 'add_animal'))
})

test_that("check_genetic() returns an error if var.ini is negative or null, or even not a number ",{
  expect_error(check_genetic(model = 'add_animal', pedigree = ped, id = id, var.ini = -1.1 ))
  expect_error(check_genetic(model = 'add_animal', pedigree = ped, id = id, var.ini = 0 ))
  expect_error(check_genetic(model = 'add_animal', pedigree = ped, id = id, var.ini = 'test' ))
})

test_that("check_genetic() returns an error if pedigree is not of class pedigree or data.frame ",{
  expect_error(check_genetic(model = 'add_animal', pedigree = FALSE, id = id, var.ini = var.ini ))
})

##Model competition##
coordinates <- matrix(c(1,2,-1,0,0,1,-1,1),4,2)
var.ini.mat <- matrix(c(1, -.5, -.5, 1), 2, 2)
pec<- list(a = FALSE, b = TRUE, c = FALSE)

test_that("The competition model runs without error",{

  ## Try alternative correct specification formats:
  comp_minimalspec <- list(
    try(check_genetic(model = 'competition',
                      pedigree = ped,   # pedigree
                      id = id,          # vector
                      coordinates = coordinates,  # matrix
                      pec = TRUE)),     # logical spec; default var.ini
    try(check_genetic(model = 'competition',
                      pedigree = ped,   # pedigree
                      id = 'id',        # variable name with data spec
                      coordinates = as.data.frame(coordinates),  # data.frame
                      pec = list(var = 1),  # list spec with var.ini abbrev.
                      var.ini = var.ini.mat,  # user var.ini
                      data = dat)),
    try(check_genetic(model = 'competition',
                      pedigree = as.data.frame(ped)[-(1:2),],  # data.frame obs.
                      id = id,          # vector
                      coordinates = as.list(as.data.frame(coordinates)),  # list
                      pec = list(present = TRUE, var.ini = 2),  # full spec
                      var.ini = var.ini.mat)), # var.ini spec
    try(check_genetic(model = 'comp',   # abbreviated name
                      pedigree = ped,   # pedigree
                      id = id,          # vector
                      coordinates = coordinates,  # matrix
                      pec = list(pres = TRUE)))  # list spec; default var.ini
  )

  for (i in seq.int(comp_minimalspec)) {
    x <- comp_minimalspec[[i]]
    var.ini.default <- c(TRUE, FALSE, FALSE, TRUE)
    
    expect_false(inherits(x, "try-error"))
    
    expect_true(setequal(names(x),
                         c('model', 'pedigree', 'id', 'coordinates', 'pec',
                           'competition_decay', 'var.ini', 'autofill')))

    expect_true(setequal(names(x$pec), c('present', 'var.ini')))
    
    ## var.ini should have been added with the default value
    ## in the cases where isTRUE(var.ini.default[[i]])
    ## and the attribute 'var.ini.default' set to TRUE
    var.ini.def <- diag(breedR.getOption('default.initial.variance'), 2)
    var.ini.def[1,2] <- var.ini.def[2,1] <- -var.ini.def[1,1]/2
    
    if (var.ini.default[[i]]) {
      expect_equal(x$var.ini, var.ini.def)
    }
    
    expect_equal(attr(x, 'var.ini.default'), var.ini.default[[i]])
  }
})

test_that("check_genetic() returns an error if missing 'coordinates' component",{
  expect_error(check_genetic(
    model = 'competition', pedigree = ped, id = id, var.ini = var.ini
  ))
})

test_that("check_genetic() returns an error if var.ini is not a SPD matrix",{
  expect_error(
    check_genetic(
      model = 'competition', pedigree = ped, coordinates = coordinates,
      id = id, var.ini = diag(8,2,3)
    )
  )
  expect_error(
    check_genetic(
      model = 'competition', pedigree = ped, coordinates = coordinates,
      id = id, var.ini = diag(-1,4,4)
    )
  )
})

test_that("check_genetic() returns an error if coordinates has not exactly two columns",{
  expect_error(
    check_genetic(
      model = 'competition', pedigree = ped, id = id, var.ini = var.ini
      , coordinates = matrix(c(1,4,6,8,5,2,3,1,5,2,1,1),4,3)
    )
  )
})

test_that("check_genetic() returns an error if pec is not a named list with logical elements",{
  expect_error(
    check_genetic(
      model = 'competition', pedigree = ped, id = id, var.ini = var.ini
      , coordinates = coordinates, pec = list(FALSE, TRUE, TRUE)
    )
  )
  expect_error(
    check_genetic(
      model = 'competition', pedigree = ped, id = id, var.ini = var.ini
      , coordinates = coordinates, pec = list(a = 5, b =
                                                'TRUE', c = TRUE)
    )
  )
})

test_that("check_genetic() returns an error if competition_decay is not a positive number",{
  expect_error(
    check_genetic(
      model = 'competition', pedigree = ped, id = id, var.ini = var.ini
      , coordinates = coordinates, pec = pec, competition_decay = -5
    )
  )
  expect_error(
    check_genetic(
      model = 'competition', pedigree = ped, id = id, var.ini = var.ini
      , coordinates = coordinates, pec = pec, competition_decay = 'test'
    )
  )
})


## Spatial


## Model splines##
var.ini <- 1.2
n.knots <- c(7,7)

test_that("Minimal specification of splines",{
  test_splines <- check_spatial(model = 'splines',
                                coordinates = coordinates)
  x <- test_splines

  expect_false(inherits(x, "try-error"))
  
  expect_true(all(names(x) %in%
                    c('model', 'coordinates', 'var.ini', 'autofill', 'sparse')))
  
  ## var.ini should have been added with the default value
  ## and the attribute 'var.ini.default' set to TRUE
  expect_equal(x$var.ini, breedR.getOption('default.initial.variance'))
  expect_true(attr(x, 'var.ini.default'))
  
})

test_that("Full specification of splines",{
  test_splines <- check_spatial(model = 'splines',
                                coordinates = coordinates,
                                n.knots = n.knots,
                                var.ini = var.ini)
  expect_false(inherits(test_splines, "try-error"))
})

test_that("check_spatial returns an error if missing 'coordinates' component",{
  expect_error(check_spatial(model = 'splines',
                             n.knots = n.knots,
                             var.ini = var.ini))
})

test_that("check_spatial returns an error if n.knots is not a vector of two integers",{
  expect_error(check_spatial(model = 'splines', coordinates = coordinates, n.knots = c(3,3,3)
                             , var.ini = var.ini))
  expect_error(check_spatial(model = 'splines', coordinates = coordinates, n.knots = TRUE
                             , var.ini = var.ini))
  expect_error(check_spatial(model = 'splines', coordinates = coordinates, n.knots = c(1.2,1.2)
                             , var.ini = var.ini))
})

## Model AR##
rho <- c(0.3,0.3)

test_that("The AR model runs without error",{
  test_ar <- check_spatial(model = 'AR', coordinates = coordinates, rho = rho, var.ini = var.ini)
  expect_false(inherits(test_ar, "try-error"))
})

test_that("check_spatial returns an error if coordinates is not a two-dimensions vector",{
  expect_error(check_spatial(model = 'AR', coordinates = diag(3), rho = rho, var.ini = var.ini))
})

test_that("check_spatial returns an error if rho does not contain what is expected",{
  expect_error(check_spatial(model = 'AR', coordinates = coordinates, rho = c(-2,1), var.ini = var.ini),
               'must contain numbers strictly between -1 and 1')
  expect_error(check_spatial(model = 'AR', coordinates = coordinates, rho = matrix(c(0.5,0,1,0),2,2), var.ini = var.ini),
               'must contain numbers strictly between -1 and 1')
  expect_error(check_spatial(model = 'AR', coordinates = coordinates, rho = c(.1,.1,.1), var.ini = var.ini),
               'must contain exactly two components')
  expect_error(check_spatial(model = 'AR', coordinates = coordinates, rho = 'test', var.ini = var.ini),
               'must be numeric')
})


## Generic

x1 <- list(inc = matrix((1:12),4,3), cov = diag(3), var.ini = 6)
x2 <- list(inc = matrix((3:8),3,2), pre = diag(2), var.ini = 4)
x <- list (a = x1, b = x2)


test_that("The function check_generic runs without error",{
  test_gen <- check_generic(x)
  expect_false(inherits(test_gen, "try-error"))
})


test_that("check_generic returns null if specification is empty",{
  expect_null(check_generic())
})

test_that("check_generic returns an error if argument x is not a list",{
  expect_error(check_generic(c(1,1)))
})

test_that("check_generic returns an error if argument x not a named list",{
  expect_error(check_generic(list(x1, x2)))
})

test_that("check_generic returns an error if all elements in x are not lists",{
  expect_error(check_generic(list(a = x1, b = x2, c = 5)))
})

test_that("check_generic returns an error if x is not a named list with different names",{
  expect_error(check_generic(list(a = x1, a = x2)))
})


test_that("check_generic returns an error if the incidence matrix is missing",{
  expect_error(check_generic(list(a = list( cov = diag(3), var.ini = 6))))
})

test_that("check_generic returns an error if both covariance and precision matrix are given",{
  expect_error(check_generic(list(a = list(inc = matrix((1:12),4,3), cov = diag(3),
                                           pre = diag(2), var.ini = 6))))
})

test_that("check_generic returns an error if both covariance and precision matrix are missing",{
  expect_error(check_generic(list(a = list(inc = matrix((1:12),4,3), var.ini = 6))))
})

test_that("check_generic returns an error if cov is not a matrix",{
  expect_error(check_generic(list(a = list(inc = matrix((1:12),4,3), cov = 'test', var.ini = 6))))
  expect_error(check_generic(list(a = list(inc = matrix((1:12),4,3), cov = c(1:3), var.ini = 6))))
})

test_that("check_generic returns an error if the dimensions are inconsistent",{
  expect_error(check_generic(list(a = list(inc = matrix((1:12),4,3), cov = diag(5), var.ini = 6))))
})

test_that("check_generic returns an error var. ini is not a positive number",{
  expect_error(check_generic(list(a = x1, b = list(inc = matrix((3:8),3,2), pre = diag(2), var.ini = -5))))
  expect_error(check_generic(list(a = x1, b = list(inc = matrix((3:8),3,2), pre = diag(2), var.ini = 'test'))))
})




context('Check Variances')

test_that('validate_variance() returns TRUE for correct variance specifications', {
  expect_true(validate_variance(1))
  expect_true(validate_variance(1.5))
  expect_true(validate_variance(1000))
  expect_true(validate_variance(matrix(c(1,-.5,-.5,1), 2, 2)))
})

test_that('validate_variance() stops for inconsistent values of variance', {
  expect_error(validate_variance(c(1, 1)), 'square matrix')
  expect_error(validate_variance(c(1, 1), dim = c(1,1)), 'square matrix')
  expect_error(validate_variance(c(1, 1), dim = c(1,2)), 'square matrix')

  expect_error(validate_variance(0), 'SPD matrix')
  expect_error(validate_variance(matrix(1, 2, 2)), 'SPD matrix')
  expect_error(validate_variance(-1.5, where = 'test'), 'SPD matrix')
  
  expect_error(validate_variance(1, dim = c(2,2)), '2x2 matrix')
})

test_that('default_initial_variance() works as expected', {
  
  ## One trait: always return half the phenotypic variance
  x <- runif(100)
  
  expect_identical(as.matrix(var(x)/2), default_initial_variance(x))
  expect_identical(as.matrix(var(x)/2), default_initial_variance(x, cor.trait = 0))
  expect_identical(as.matrix(var(x)/2), default_initial_variance(x, cor.effect = 0))
  
  ## One trait - 2 dimensional effect (e.g. competition)
  x <- runif(100)
  dim = 2
  default.covar <- 0.1*var(x)/2
  
  div <- default_initial_variance(x, dim = dim)
  expect_identical(rep(as.matrix(var(x)/2), dim), diag(div))
  expect_identical(default.covar, div[2,1])
  expect_identical(default.covar, div[2,1])
  
  div <- default_initial_variance(x, dim = dim, cor.effect = 0)
  expect_identical(rep(as.matrix(var(x)/2), dim), diag(div))
  expect_identical(0, div[2,1])
  
  ## 5 traits: half phenotypic variance of each trait and default covariances
  ## unless diag
  x <- matrix(runif(500), ncol=5)
  default.covar <- 0.1*gmean(diag(var(x)))/2
  
  div <- default_initial_variance(x)
  expect_identical(diag(var(x))/2, diag(div))
  expect_equal(default.covar, div[2,1])
  
  divd <- default_initial_variance(x, cor.trait = 0)
  expect_identical(diag(var(x))/2, diag(divd))
  expect_identical(0, divd[2, 1])
  
  ## Fail at constant traits
  x[, 5] <- 123
  expect_error(default_initial_variance(x), 'Trait 5 is constant.')
})
