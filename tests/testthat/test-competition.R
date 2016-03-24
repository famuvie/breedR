old.op <- options(warn = -1,  # suppressWarnings
                  show.error.messages = FALSE)  # silent try
on.exit(options(old.op))


#### Context: competition infrastructure ####
context("competition infrastructure")

## Minimal dataset
dat <- data.frame(id   = 1:6,
                  sire = c(11, 11, 2, 3, 11, 3),
                  dam  = c(12, NA, 1, 12, 12, 1),
                  x    = c(1,2,-1,0,0,1),
                  y    = c(-1,0,0,1,-1,1))
## Corresponding pedigree with additional offspring
ped <- build_pedigree(1:3, data = rbind(dat, c(7, 1, 2)))
var.ini.mat <- matrix(c(1, -.5, -.5, 1), 2, 2)


test_that("Valid alternative model specifications pass check_genetic()", {
  ## specify var.ini, but use pec=FALSE (as by default)
  expect_error(
    check_genetic(model = 'competition',
                  pedigree = ped,
                  coordinates = dat[, c('x', 'y')],
                  id = dat$id,
                  var.ini = var.ini.mat,
                  response = rnorm(nrow(dat))),
    NA
  )
  
  ## non-recoded pedigree
  idx <- attr(ped, 'map')[dat$id]
  dat[, 1:3] <- as.data.frame(ped)[idx, ]
  expect_error(
    check_genetic(model = 'competition',
                  pedigree = dat[, 1:3],
                  coordinates = dat[, c('x', 'y')],
                  id = dat$id,
                  var.ini = var.ini.mat,
                  response = rnorm(nrow(dat))),
    NA
  )
  
})



test_that("Invalid alternative model specifications fail check_genetic()", {

  ## incomplete non-recoded pedigrees
  idx <- attr(ped, 'map')[dat$id]
  dat[, 1:3] <- as.data.frame(ped)[idx, ]
  expect_error(
    check_genetic(model = 'competition',
                  pedigree = dat[-nrow(dat), 1:3],
                  coordinates = dat[, c('x', 'y')],
                  id = dat$id,
                  var.ini = var.ini.mat,
                  response = rnorm(nrow(dat))), 
    'The following individuals in id are not represented'
  )
})



test_that("additive_genetic_competition() works as expected", {
  ## Full specification, from minimal input
  comp.spec <- check_genetic(
    model = 'competition',
    pedigree = ped,
    id = dat$id,
    coordinates = dat[, c('x', 'y')],
    pec = TRUE,
    response = rnorm(nrow(dat))
  )
  
  res <- with(
    comp.spec,
    additive_genetic_competition(
      pedigree    = pedigree,
      coordinates = coordinates,
      id          = id,
      decay       = competition_decay,
      autofill    = autofill
    )
  )
  
  expect_is(res, c("additive_genetic_competition", "additive_genetic", "genetic", 
                   "competition", "spatial", "random", "breedr_effect"))
  expect_equal(length(res), 5)
  # Incidence matrix
  inc.mat <- model.matrix(res)
  expect_is(inc.mat, 'sparseMatrix') # a permutation Matrix
  expect_equal(nrow(inc.mat), nrow(dat))
  # Covariance matrix
  cov.mat <- get_structure(res)
  expect_is(cov.mat, 'sparseMatrix') 
  expect_equal(ncol(cov.mat), nrow(as.data.frame(ped)))
  expect_equal(ncol(inc.mat), nrow(cov.mat))
})
