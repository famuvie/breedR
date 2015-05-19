old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

### Tests of functions for checking model components ###

context("check_genetic")

##Model add_animal##
dat <- data.frame(id = 1:4,
                  sire = c(11, 11, 2, 3),
                  dam  = c(12, NA, 1, 12))
ped <- build_pedigree(1:3, data = dat)
id <- dat$id
var.ini <- 1.5

test_that("The add_animal model runs without error",{
  test_add_animal <- check_genetic(model = 'add_animal', pedigree = ped, id = id, var.ini = var.ini )
  expect_that(!inherits(test_add_animal, "try-error"), is_true())
})

test_that("check_genetic returns an error if missing values",{
  expect_error(check_genetic(pedigree = ped, id = id, var.ini = var.ini ))
  expect_error(check_genetic(model = 'add_animal', pedigree = ped, var.ini = var.ini ))
  expect_error(check_genetic(model = 'add_animal'))
})

test_that("check_genetic returns an error if var.ini is negative or null, or even not a number ",{
  expect_error(check_genetic(model = 'add_animal', pedigree = ped, id = id, var.ini = -1.1 ))
  expect_error(check_genetic(model = 'add_animal', pedigree = ped, id = id, var.ini = 0 ))
  expect_error(check_genetic(model = 'add_animal', pedigree = ped, id = id, var.ini = 'test' ))
})

test_that("check_genetic returns an error if pedigree is not of class pedigree or data.frame ",{
  expect_error(check_genetic(model = 'add_animal', pedigree = FALSE, id = id, var.ini = var.ini ))
})

##Model competition##
coordinates <- matrix(c(1,2,-1,0,0,1,-1,1),4,2)
var.ini <- diag(4)
pec<- list(a = FALSE, b = TRUE, c = FALSE)

test_that("The competition model runs without error",{
  test_compet <- check_genetic(model = 'competition', pedigree = ped, id = id, 
                               coordinates = coordinates, pec= pec, var.ini = var.ini )
  expect_that(!inherits(test_compet, "try-error"), is_true())
})

test_that("check_genetic returns an error if missing 'coordinates' component",{
  expect_error(check_genetic(model = 'competition', pedigree = ped, id=id, var.ini = var.ini ))
})

test_that("check_genetic returns an error if var.ini is not a SPD matrix",{
  expect_error(check_genetic(model = 'competition', pedigree = ped, coordinates = coordinates,
                             id=id, var.ini = diag(8,2,3)))
  expect_error(check_genetic(model = 'competition', pedigree = ped, coordinates = coordinates,
                             id=id, var.ini = diag(-1,4,4)))
})

test_that("check_genetic returns an error if coordinates has not exactly two columns",{
  expect_error(check_genetic(model = 'competition', pedigree = ped, id=id, var.ini = var.ini
                             , coordinates = matrix(c(1,4,6,8,5,2,3,1,5,2,1,1),4,3)))
})

test_that("check_genetic returns an error if pec is not a named list with logical elements",{
  expect_error(check_genetic(model = 'competition', pedigree = ped, id=id, var.ini = var.ini
                             , coordinates = coordinates, pec = list(FALSE, TRUE, TRUE)))
  expect_error(check_genetic(model = 'competition', pedigree = ped, id=id, var.ini = var.ini
                             , coordinates = coordinates, pec = list(a= 5, b='TRUE', c=TRUE)))
})

test_that("check_genetic returns an error if competition_decay is not a positive number",{
  expect_error(check_genetic(model = 'competition', pedigree = ped, id=id, var.ini = var.ini
                             , coordinates = coordinates, pec = pec, competition_decay = -5))
  expect_error(check_genetic(model = 'competition', pedigree = ped, id=id, var.ini = var.ini
                             , coordinates = coordinates, pec = pec, competition_decay = 'test'))
})



context("check_spatial")

##Model splines##
var.ini <- 1.2
n.knots <- c(7,7)

test_that("The splines model runs without error",{
  test_splines <- check_spatial(model = 'splines', coordinates = coordinates, n.knots = n.knots
                                , var.ini = var.ini)
  expect_that(!inherits(test_splines, "try-error"), is_true())
})

test_that("check_spatial returns an error if missing 'coordinates' component",{
  expect_error(check_spatial(model = 'splines', n.knots = n.knots
                             , var.ini = var.ini))
})

test_that("check_spatial returns an error if var.ini is not a positive number",{
  expect_error(check_spatial(model = 'splines', coordinates = coordinates, n.knots = n.knots
                             , var.ini = -1))
  expect_error(check_spatial(model = 'splines', coordinates = coordinates, n.knots = n.knots
                             , var.ini = c(3,3)))
})

test_that("check_spatial returns an error if n.knots is not a vector of two integers",{
  expect_error(check_spatial(model = 'splines', coordinates = coordinates, n.knots = c(3,3,3)
                             , var.ini = var.ini))
  expect_error(check_spatial(model = 'splines', coordinates = coordinates, n.knots = TRUE
                             , var.ini = var.ini))
  expect_error(check_spatial(model = 'splines', coordinates = coordinates, n.knots = c(1.2,1.2)
                             , var.ini = var.ini))
})

##Model AR##
rho <- c(0.3,0.3)

test_that("The AR model runs without error",{
  test_ar <- check_spatial(model = 'AR', coordinates = coordinates, rho = rho, var.ini = var.ini)
  expect_that(!inherits(test_ar, "try-error"), is_true())
})

test_that("check_spatial returns an error if missing 'var.ini' component",{
  expect_error(check_spatial(model = 'AR', coordinates = coordinates, rho = rho))
})

test_that("check_spatial returns an error if coordinates is not a two-dimensions vector",{
  expect_error(check_spatial(model = 'AR', coordinates = diag(3), rho = rho, var.ini = var.ini))
})

test_that("check_spatial returns an error if rho does not contain what is expected",{
  expect_error(check_spatial(model = 'AR', coordinates = coordinates, rho = c(1,1,1), var.ini = var.ini))
  expect_error(check_spatial(model = 'AR', coordinates = coordinates, rho = c(-2,1), var.ini = var.ini))
  expect_error(check_spatial(model = 'AR', coordinates = coordinates, rho = diag(0.5,2,2), var.ini = var.ini))
  expect_error(check_spatial(model = 'AR', coordinates = coordinates, rho = 'test', var.ini = var.ini))
})

context("check_generic")

x1 <- list(inc = matrix((1:12),4,3), cov = diag(3), var.ini = 6)
x2 <- list(inc = matrix((3:8),3,2), pre = diag(2), var.ini = 4)
x <- list (a = x1, b = x2)

test_gen <- check_generic(x)

test_that("The function check_generic runs without error",{
  expect_that(!inherits(test_gen, "try-error"), is_true())
})


test_that("check_generic returns an error if argument x is missing",{
  expect_error(check_generic())
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

## TODO : fix what causes an error in the following instructions 

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
