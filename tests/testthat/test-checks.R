old.op <- options(warn = -1,  # suppressWarnings
                  show.error.messages = FALSE)  # silent try
on.exit(options(old.op))

### Tests of functions for checking model components ###

context("Checks for model components")

## Model add_animal ##
dat <- data.frame(id   = 1:4,
                  sire = c(11, 11, 2, 3),
                  dam  = c(12, NA, 1, 12),
                  y    = rnorm(4),
                  z    = rnorm(4, sd = 2))
ped <- build_pedigree(1:3, data = dat)
id <- dat$id
var.ini <- 1.5

divf <- breedR.getOption('default.initial.variance')

test_that("The minimal add_animal specification checks without error and completes missing values",{
  
  ## Try alternative correct specification formats:
  animal_minimalspec <- list(
    list(model = 'add_animal',
         pedigree = ped,
         id = id),
    list(model = 'add_animal',
         pedigree = ped,
         id = 'id',
         data = dat),
    list(model = 'add_animal',
         pedigree = as.data.frame(ped)[-(1:2),],
         id = id),
    list(model = 'add',
         pedigree = ped,
         id = id)
  )

  ## Single-trait
  animal_check_input <- lapply(animal_minimalspec, c, response = list(dat$y))
  animal_checks <- 
    lapply(animal_check_input, function(x) try(do.call('check_genetic', x)))
  
  for (x in animal_checks) {
    expect_false(inherits(x, "try-error"))
    
    expect_true(setequal(names(x),
                         c('model', 'pedigree', 'id', 'var.ini', 'autofill')))
    
    ## var.ini should have been added with the default value
    ## and the attribute 'var.ini.default' set to TRUE
    expect_equal(x$var.ini, eval(divf)(dat$y, digits = 2))
    expect_true(attr(x, 'var.ini.default'))
  }

  ## Two-trait
  animal_check_input <- lapply(animal_minimalspec, c, 
                               response = list(cbind(dat$y, dat$z)))
  animal_checks <- 
    lapply(animal_check_input, function(x) try(do.call('check_genetic', x)))
  
  for (x in animal_checks) {
    expect_false(inherits(x, "try-error"))
    
    expect_true(setequal(names(x),
                         c('model', 'pedigree', 'id', 'var.ini', 'autofill')))
    
    ## var.ini should have been added with the default value
    ## and the attribute 'var.ini.default' set to TRUE
    expect_equal(x$var.ini, eval(divf)(dat[, c('y', 'z')], digits = 2))
    expect_true(attr(x, 'var.ini.default'))
  }
})

test_that("check_genetic() returns an error if missing arguments",{
  expect_error(check_genetic(pedigree = ped, id = id, var.ini = var.ini),
               'model required')
  expect_error(check_genetic(model = 'add_animal', pedigree = ped, var.ini = var.ini),
               'id required')
  expect_error(check_genetic(model = 'add_animal'), 'pedigree required')
})

test_that("check_genetic() returns an error if var.ini is inconsistent",{
  check_var <- function(x, Y = dat$y)
    check_genetic(model = 'add_animal', 
                  pedigree = ped, 
                  id = id, 
                  var.ini = x, 
                  response = Y)

  ## var.ini not of right dimension
  expect_error(check_var(1, Y = dat[, c('y', 'z')]), '2x2 matrix')
})

test_that("check_genetic() returns an error if pedigree is not of class pedigree or data.frame ",{
  expect_error(check_genetic(model = 'add_animal', 
                             pedigree = FALSE, 
                             id = id, 
                             var.ini = var.ini),
               'argument pedigree')
})

## Model competition ##
coordinates <- matrix(c(1,2,-1,0,0,1,-1,1),4,2)
var.ini.mat <- matrix(c(1, -.5, -.5, 1), 2, 2)
pec<- list(a = FALSE, b = TRUE, c = FALSE)

test_that("Correctly-specified competition models runs without error",{

  ## Try alternative correct specification formats:
  comp_minimalspec <- list(
    list(model = 'competition',
         pedigree = ped,   # pedigree
         id = id,          # vector
         coordinates = coordinates,  # matrix
         pec = TRUE),     # logical spec; default var.ini
    list(model = 'competition',
         pedigree = ped,   # pedigree
         id = 'id',        # variable name with data spec
         coordinates = as.data.frame(coordinates),  # data.frame
         pec = list(var = 1),  # list spec with var.ini abbrev.
         var.ini = var.ini.mat,  # user var.ini
         data = dat),
    list(model = 'competition',
         pedigree = as.data.frame(ped)[-(1:2),],  # data.frame obs.
         id = id,          # vector
         coordinates = as.list(as.data.frame(coordinates)),  # list
         pec = list(present = TRUE, var.ini = 2),  # full spec
         var.ini = var.ini.mat), # var.ini spec
    list(model = 'comp',   # abbreviated name
         pedigree = ped,   # pedigree
         id = id,          # vector
         coordinates = coordinates,  # matrix
         pec = list(pres = TRUE))  # list spec; default var.ini
  )

  ## Single-trait
  comp_check_input <- lapply(comp_minimalspec, c, response = list(dat$y))

  comp_checks <- 
    lapply(comp_check_input, function(x) do.call('check_genetic', x))
  
  all.names <- c('model', 'pedigree', 'id', 'coordinates', 'pec',
                 'competition_decay', 'var.ini', 'autofill')
  expect_true(all(sapply(lapply(comp_checks, names), setequal, all.names)))
  
  pec.names <- lapply(lapply(comp_checks, function(x) x$pec), names)
  expect_true(all(sapply(pec.names, setequal, c('present', 'var.ini'))))

  ## var.ini should have been added with the default value
  ## in the cases where isTRUE(var.ini.default[[i]])
  ## and the attribute 'var.ini.default' set to TRUE
  expect_defvar <- eval(divf)(dat$y, dim = 2, digits = 2)
  expect_var <- list(expect_defvar, var.ini.mat, var.ini.mat, expect_defvar)
  expect_identical(lapply(comp_checks, function(x) x$var.ini),
                   expect_var)
  
  expect_var.ini.default <- c(TRUE, FALSE, FALSE, TRUE)
  expect_identical(sapply(comp_checks, attr, 'var.ini.default'),
                   expect_var.ini.default)
  

  ## Two traits
  var.ini.mat <- Matrix::bdiag(list(var.ini.mat, var.ini.mat))
  comp_minimalspec[[2]]$var.ini <- comp_minimalspec[[3]]$var.ini <- var.ini.mat
  comp_check_input <- lapply(comp_minimalspec, c,
                             response = list(dat[, c('y', 'z')]))
  comp_checks <-
    lapply(comp_check_input, function(x) do.call('check_genetic', x))
  
  all.names <- c('model', 'pedigree', 'id', 'coordinates', 'pec',
                 'competition_decay', 'var.ini', 'autofill')
  expect_true(all(sapply(lapply(comp_checks, names), setequal, all.names)))
  
  pec.names <- lapply(lapply(comp_checks, function(x) x$pec), names)
  expect_true(all(sapply(pec.names, setequal, c('present', 'var.ini'))))
  
  ## var.ini should have been added with the default value
  ## in the cases where isTRUE(var.ini.default[[i]])
  ## and the attribute 'var.ini.default' set to TRUE
  expect_defvar <- eval(divf)(dat[, c('y', 'z')], dim = 2, digits = 2)
  expect_var <- list(expect_defvar, var.ini.mat, var.ini.mat, expect_defvar)
  expect_identical(lapply(comp_checks, function(x) x$var.ini),
                   expect_var)
  
  expect_var.ini.default <- c(TRUE, FALSE, FALSE, TRUE)
  expect_identical(sapply(comp_checks, attr, 'var.ini.default'),
                   expect_var.ini.default)
  
})

test_that("check_genetic() returns an error if missing arguments",{
  expect_error(check_genetic(model = 'competition', 
                             pedigree = ped, 
                             id = id, 
                             var.ini = var.ini.mat, 
                             response = dat$y),
               'coordinates required')
})

test_that("check_genetic() returns an error if var.ini is incorrect",{
  
  ## Single trait
  expect_error(
    check_genetic(
      model = 'competition', pedigree = ped, coordinates = coordinates,
      id = id, var.ini = diag(-1,4,4), response = dat$y
    ),
    '2x2 matrix'
  )

  ## Two traits
  expect_error(
    check_genetic(
      model = 'competition', pedigree = ped, coordinates = coordinates,
      id = id, var.ini = diag(1,2,2), response = dat[, c('y', 'z')]
    ),
    '4x4 matrix'
  )
})

test_that("check_genetic() returns an error if coordinates has not exactly two columns",{
  expect_error(
    check_genetic(
      model = 'competition', pedigree = ped, id = id, var.ini = var.ini.mat
      , coordinates = matrix(c(1,4,6,8,5,2,3,1,5,2,1,1),4,3), response = dat$y
    ),
    'two dimensions admitted for coordinates'
  )
})

test_that("check_genetic() returns an error if pec is not a named list with logical elements",{
  expect_error(
    check_genetic(
      model = 'competition', 
      pedigree = ped, 
      id = id, 
      var.ini = var.ini.mat,
      coordinates = coordinates, 
      pec = list(FALSE, TRUE, TRUE),
      response = dat$y
    ),
    'pec must be a named list'
  )
  
  expect_error(
    check_genetic(
      model = 'competition', 
      pedigree = ped, 
      id = id, 
      var.ini = var.ini.mat ,
      coordinates = coordinates, 
      pec = list(a = 5, b = 'TRUE', c = TRUE), 
      response = dat$y
    ),
    'should be one of'
  )
})

test_that("check_genetic() returns an error if competition_decay is not a positive number",{
  expect_error(
    check_genetic(
      model = 'competition', pedigree = ped, id = id, var.ini = var.ini.mat, 
      coordinates = coordinates, pec = 1, competition_decay = -5, response = dat$y
    ),
    'competition_decay > 0'
  )
  
  expect_error(
    check_genetic(
      model = 'competition', pedigree = ped, id = id, var.ini = var.ini.mat, 
      coordinates = coordinates, pec = 1, competition_decay = 'test', response = dat$y
    ),
    'is.numeric\\(competition_decay\\)'
  )
})


## Spatial
var.ini <- 1.2


## Model splines ##
n.knots <- c(7,7)


test_that("Minimal correct specification of splines",{
  
  ## One trait
  spl_check <- check_spatial(model = 'splines',
                             coordinates = coordinates,
                             response = dat$y)
  
  expect_false(inherits(spl_check, "try-error"))
  
  expect_true(all(names(spl_check) %in%
                    c('model', 'coordinates', 'var.ini', 'autofill', 'sparse')))
  
  ## var.ini should have been added with the default value
  ## and the attribute 'var.ini.default' set to TRUE
  expect_equal(spl_check$var.ini, eval(divf)(dat$y, digits = 2))
  expect_true(attr(spl_check, 'var.ini.default'))
  
  ## Two traits
  spl_check <- check_spatial(model = 'splines',
                             coordinates = coordinates,
                             response = dat[, c('y', 'z')])
  
  expect_false(inherits(spl_check, "try-error"))
  
  expect_true(all(names(spl_check) %in%
                    c('model', 'coordinates', 'var.ini', 'autofill', 'sparse')))
  
  ## var.ini should have been added with the default value
  ## and the attribute 'var.ini.default' set to TRUE
  expect_equal(spl_check$var.ini, eval(divf)(dat[, c('y', 'z')], digits = 2))
  expect_true(attr(spl_check, 'var.ini.default'))
  
})

test_that("Full correct specification of splines",{
  
  ## One trait
  spl_check <- check_spatial(model = 'splines',
                             coordinates = coordinates,
                             n.knots = n.knots,
                             var.ini = var.ini,
                             response = dat$y)
  expect_false(inherits(spl_check, "try-error"))

  ## Two traits
  spl_check <- check_spatial(model = 'splines',
                             coordinates = coordinates,
                             n.knots = n.knots,
                             var.ini = diag(rep(var.ini, 2)),
                             response = dat[, c('y', 'z')])
  expect_false(inherits(spl_check, "try-error"))
})

test_that("check_spatial() errors if 'coordinates' is wrongly specified",{
  
  expect_error(check_spatial(model = 'splines',
                             n.knots = n.knots,
                             var.ini = var.ini,
                             response = dat$y),
               'coordinates required')

  expect_error(check_spatial(model = 'splines',
                             coordinates = diag(3),
                             n.knots = n.knots,
                             var.ini = var.ini,
                             response = dat$y),
               'Only two dimensions admitted for coordinates')
})

test_that("check_spatial() errors if 'n.knots' is wrongly specified",{

  expect_nk_error <- function(x) 
    eval(bquote(
      expect_error(
        check_spatial(model = 'splines', 
                      coordinates = coordinates,
                      n.knots = .(x),
                      var.ini = var.ini,
                      response = dat$y), 
        'n.knots must be a vector of two integers')
    ))
  
  expect_nk_error(c(3,3,3))
  expect_nk_error(TRUE)
  expect_nk_error(c(1.2,1.2))
  
})

test_that("check_spatial() errors if var.ini is inconsistent",{
  
  ## Single trait
  expect_error(
    check_spatial(
      model = 'splines', coordinates = coordinates,
      id = id, var.ini = diag(-1,4,4), response = dat$y
    ),
    '1x1 matrix'
  )

  ## Two traits
  expect_error(
    check_spatial(
      model = 'splines', coordinates = coordinates,
      id = id, var.ini = diag(1,4,4), response = dat[, c('y', 'z')]
    ),
    '2x2 matrix'
  )
})


## Model AR ##
rho <- c(0.3,0.3)

test_that("The AR model runs without error",{
  
  ## One trait
  ar_check <- check_spatial(model = 'AR',
                           coordinates = coordinates, 
                           rho = rho, 
                           response = dat$y)
  
  expect_false(inherits(ar_check, "try-error"))
  expect_true(all(names(ar_check) %in%
                    c('model', 'coordinates', 'rho', 'var.ini', 'autofill', 'sparse')))
  
  ## var.ini should have been added with the default value
  ## and the attribute 'var.ini.default' set to TRUE
  expect_equal(ar_check$var.ini, eval(divf)(dat$y, digits = 2))
  expect_true(attr(ar_check, 'var.ini.default'))
  
  ## Two traits
  ar_check <- check_spatial(model = 'AR',
                            coordinates = coordinates, 
                            rho = rho, 
                            response = dat[, c('y', 'z')])
  
  expect_false(inherits(ar_check, "try-error"))
  expect_true(all(names(ar_check) %in%
                    c('model', 'coordinates', 'rho', 'var.ini', 'autofill', 'sparse')))
  
  ## var.ini should have been added with the default value
  ## and the attribute 'var.ini.default' set to TRUE
  expect_equal(ar_check$var.ini, eval(divf)(dat[, c('y', 'z')], digits = 2))
  expect_true(attr(ar_check, 'var.ini.default'))
  
})

test_that("check_spatial() errors if 'coordinates' is wrongly specified",{

  expect_error(check_spatial(model = 'AR', 
                             rho = rho, 
                             var.ini = var.ini),
               'coordinates required')

  expect_error(check_spatial(model = 'AR', 
                             coordinates = diag(3),
                             rho = rho, 
                             var.ini = var.ini),
               'Only two dimensions admitted for coordinates')
})

test_that("check_spatial() errors if 'rho' is incorrectly specified",{
  
  expect_rho <- function(x, msg) {
    eval(bquote(
      expect_error(
        check_spatial(model = 'AR', 
                      coordinates = coordinates, 
                      rho = .(x), 
                      var.ini = var.ini, 
                      response = dat$y),
        msg)
    ))
  }
  
  expect_rho(c(-2,1), 'strictly between -1 and 1')
  expect_rho(matrix(c(0.5,0,1,0),2,2), 'strictly between -1 and 1')
  expect_rho(c(.1,.1,.1), 'exactly two components')
  expect_rho('test', 'be numeric')

})

test_that("check_spatial() errors if var.ini is inconsistent",{
  
  ## Single trait
  expect_error(
    check_spatial(
      model = 'AR', coordinates = coordinates,
      id = id, var.ini = diag(-1,4,4), response = dat$y
    ),
    '1x1 matrix'
  )
  
  ## Two traits
  expect_error(
    check_spatial(
      model = 'AR', coordinates = coordinates,
      id = id, var.ini = diag(1,4,4), response = dat[, c('y', 'z')]
    ),
    '2x2 matrix'
  )
})



## Generic

x1 <- list(inc = matrix((1:12),4,3), cov = diag(3), var.ini = 6)
x2 <- list(inc = matrix((3:8),3,2), pre = diag(2), var.ini = 4)
x <- list (a = x1, b = x2)


test_that("Correct specification of individual generic elements", {

  expect_ok <- function(x, Y) {
    eval(bquote(
      expect_error(do.call('validate_generic_element',
                           c(.(x), response = list(.(Y)))),
                   NA)
    ))
  }
  
  ## Single trait
  expect_ok(x1, dat$y)
  expect_ok(x2, dat$y)
  expect_ok(x1[-3], dat$y)
  expect_ok(x2[-3], dat$y)
  
  ## Two traits
  x1$var.ini <- x2$var.ini <- diag(1,2)
  expect_ok(x1, dat[, c('y', 'z')])
  expect_ok(x2, dat[, c('y', 'z')])
  expect_ok(x1[-3], dat[, c('y', 'z')])
  expect_ok(x2[-3], dat[, c('y', 'z')])
  
})

test_that("Correct specifications of check_generic()",{
  expect_error(check_generic(x, response = dat$y), NA)
})


test_that("check_generic returns null if specification is empty",{
  expect_null(check_generic())
})

test_that("check_generic() errors if specification is wrong",{
  
  expect_error(check_generic(c(1,1)), 'be a list')
  
  expect_error(check_generic(list(x1, x2)), 'be a named list')
  
  expect_error(check_generic(list(a = x1, b = x2, c = 5)), 'be list elements')
  
  expect_error(check_generic(list(a = x1, a = x2)), 'different names')
  
  expect_error(check_generic(list(a = list(cov = diag(3), var.ini = 6)), 
                             response = dat$y), 'incidence required')
  
  expect_error(check_generic(list(a = list(inc = matrix((1:12),4,3), 
                                           cov = diag(3),
                                           pre = diag(2))), 
                             response = dat$y), 
               'one argument between covariance and precision')
  
  expect_error(check_generic(list(a = list(inc = matrix((1:12),4,3), var.ini = 6)), 
                             response = dat$y), 
               'one argument between covariance and precision')

  expect_error(check_generic(list(a = list(inc = matrix((1:12),4,3), cov = 'test')), 
                             response = dat$y),
               'covariance must be of type matrix')
  
  expect_error(check_generic(list(a = list(inc = matrix((1:12),4,3), cov = c(1:3))), 
                             response = dat$y),
               'covariance must be of type matrix')

  expect_error(check_generic(list(a = list(inc = matrix((1:12),4,3), cov = diag(5))), 
                             response = dat$y),
               'conformant incidence and covariance')
})




context('Check Variances')

test_that('validate_variance() returns TRUE for correct variance specifications', {
  expect_true(validate_variance(1))
  expect_true(validate_variance(1.5))
  expect_true(validate_variance(1000))
  expect_true(validate_variance(matrix(c(1,-.5,-.5,1), 2, 2)))
  
  ## sparse matrices
  expect_true(validate_variance(Matrix::bdiag(matrix(c(1,-.5,-.5,1), 2, 2))))
  
  ## named (possibly differently) matrices
  testM <- structure(matrix(c(1,-.5,-.5,1),2,2),
                     dimnames = list(1:2, 3:4))
  expect_true(validate_variance(testM))
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
  expect_equal(rep(as.matrix(var(x)/2), dim), diag(div))
  expect_equal(default.covar, div[2,1])
  expect_equal(default.covar, div[1,2])
  
  div <- default_initial_variance(x, dim = dim, cor.effect = 0)
  expect_identical(rep(as.matrix(var(x)/2), dim), diag(div))
  expect_identical(0, div[2,1])
  
  ## 5 traits: half phenotypic variance of each trait and default covariances
  ## unless diag
  x <- matrix(runif(500), ncol=5)
  
  expect_equal(default_initial_variance(x), var(x)/2)

  expect_equal(default_initial_variance(x, cor.trait = 0),
               diag(diag(var(x)/2)))

  ## Fail at constant traits
  x[, 5] <- 123
  expect_error(default_initial_variance(x), 'Trait 5 is constant.')
  
  ## 2 traits - 2 dimensions: positive definiteness
  x <- data.frame(y=rnorm(4), z=rnorm(4, sd = 2))
  varx <- default_initial_variance(x, dim = 2)
  expect_error(validate_variance(varx), NA)
  
})

