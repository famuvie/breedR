suppressPackageStartupMessages(
 require(spam)
)  

#### Context: Models with several effects working together ####
context("Models with several effects working together") 

# TODO:
# mixed effects models * spatial models * genetic effects

dat <- try(
  breedR.sample.phenotype(
    fixed   = c(mu = 10, x = 2),
    random = list(u = list(nlevels = 3,
                           sigma2  = 1)),
    genetic = list(model    = 'competition',
                   Nparents = c(10, 10),
                   sigma2_a = matrix(c(2, -1, -1, 2), 2, 2),
                   competition_decay = 1,
                   check.factorial = FALSE,
                   pec = 0.5),
    spatial = list(model     = 'splines',
                   grid.size = c(10, 10),
                   n.knots   = c(3, 3),
                   sigma2_s  = 1),
    residual.variance = 1)[-(1:20), ]
)

test_that("Simulation of splines+competition succeeds", {
  expect_false(inherits(dat, 'try-error'))
})

if (!inherits(dat, 'try-error')) {
  coord <- dat[, c('Var1', 'Var2')]
  
  res <- try(
    suppressMessages(
      remlf90(
        phenotype ~ X.x,
        random = ~ u,
        genetic = list(model = 'competition',
                       pedigree = dat[, 1:3],
                       id = 'self',
                       coord = coord,
                       competition_decay = 1,
                       pec = list(present = TRUE)),
        spatial = list(model = 'splines',
                       coord = coord,
                       n.knots = c(3, 3)),
        method = 'em',  # ai seems to work in interactive mode, but fails to converge in test()
        data = dat)
    )
  )
  
  test_that("Fitting of splines+competition simulation succeeds", {
    expect_false(inherits(res, 'try-error'))
  })
  
  if (!inherits(res, 'try-error')) {
    
    test_that("coordinates.breedR() extracts and unifies multiple equal sets of coordinates", {
      expect_equal(coordinates(res), coord, check.attributes = FALSE)
      
      ## Introduce a difference in the coordinates sets and verify that it is detected
      res_fail <- res
      res_fail$effects$spatial$effects[[1]]$coordinates[1,1] <- 0
      expect_error(coordinates(res_fail), 'different coordinate sets in the model')
    })
    
  }
}
