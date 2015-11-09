old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))


#### Context: additive_genetic_competition() ####
context("competition infrastructure")


#' dat <- data.frame(id   = 1:5,
#'                   sire = c(11, 11, 2, 3, 2),
#'                   dam  = c(12, NA, 1, 12, 1),
#'                   x    = c(rep(1:2, times = 2), 3),
#'                   y    = c(rep(1:2, each = 2), 3))
#' ped <- build_pedigree(1:3, data = dat)
#' breedR:::additive_genetic_competition(ped, coord = dat[, c('x', 'y')], dat$id, 2)


## Minimal dataset
dat <- data.frame(id   = 1:6,
                  sire = c(11, 11, 2, 3, 11, 3),
                  dam  = c(12, NA, 1, 12, 12, 1),
                  x    = c(1,2,-1,0,0,1),
                  y    = c(-1,0,0,1,-1,1))
ped <- build_pedigree(1:3, data = rbind(dat, c(7, 1, 2)))
var.ini <- 1.5
var.ini.mat <- matrix(c(1, -.5, -.5, 1), 2, 2)

## Full specification, from minimal input
comp.spec <- check_genetic(model = 'competition',
                           pedigree = ped,
                           id = dat$id,
                           coordinates = dat[, c('x', 'y')],
                           pec = TRUE)

res <- try(with(comp.spec,
                additive_genetic_competition(pedigree    = pedigree,
                                             coordinates = coordinates,
                                             id          = id,
                                             decay       = competition_decay,
                                             autofill    = autofill)))

test_that("additive_genetic_competition() works as expected", {
  expect_false(inherits(res, 'try-error'))
  
  if (!inherits(res, 'try-error')) {
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
  }
})
