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
test_that("The competition model runs without error",{
  test_compet <- check_genetic(model = 'competition', pedigree = ped, id = id, 
                               coordinates = coordinates, var.ini = var.ini )
  expect_that(!inherits(test_compet, "try-error"), is_true())
})

