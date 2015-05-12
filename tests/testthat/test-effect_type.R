## Test of the function effect_type 

old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))


context("effect_type")

test.inc8x4 <- as(rep(1:4, 2), 'indMatrix')   # An 8x4 incidence matrix
test.cov4x4 <- with(L <- Matrix::tril(matrix(sample(16),4)),
                    Matrix::t(L)%*%L)         # a SPD 4x4 covariance matrix

##Test with an object of type generic
gtest <- generic(incidence  = test.inc8x4, covariance = test.cov4x4)
test_that("The function returns the type expected",{
  expect_output(effect_type(gtest),"random")
})

##Test with an object of type splines
x.loc <- 1:100
y.loc <- seq(1000, by = 5, length = 51)
coord <- expand.grid(x.loc, y.loc)
result <- splines(coord)
test_that("The function returns the type expected",{
  expect_output(effect_type(result),"random")
})

##Test with an object of type competition
x = c(rep(1:2, times = 2), 3)
y = c(rep(1:2, each = 2), 3)
dat <- data.frame(id   = 1:5,
                  sire = c(11, 11, 2, 3, 2),
                  dam  = c(12, NA, 1, 12, 1),
                  x    = x,
                  y    = y)

coordinates <- dat[, c('x', 'y')]
covariance <- diag(nrow(coordinates))
comptest <- competition(coordinates,covariance, decay=1,autofill=TRUE)
test_that("The function returns the type expected",{
  expect_output(effect_type(comptest),"random")
})

##Test with an object of type additive_genetic_competition
ped <- build_pedigree(1:3, data = dat)
addtest <- additive_genetic_competition(ped, coord = coordinates, dat$id, 1, autofill=TRUE)
test_that("The function returns the type expected",{
  expect_output(effect_type(addtest),"random")
})

##Test with an object of type additive_genetic
ped2 <- pedigreemm::pedigree(sire = c(NA,NA,1, 1,4,5),
                            dam  = c(NA,NA,2,NA,3,2),
                            label= 1:6)
inc <- cbind(0, 1, diag(4))
addtest2 <- additive_genetic(ped2, inc)
test_that("The function returns the type expected",{
  expect_output(effect_type(addtest2),"random")
})

##Test with an object of type additive_genetic_animal
dat <- data.frame(id = 1:4,
                  sire = c(11, 11, 2, 3),
                  dam  = c(12, NA, 1, 12))
ped3 <- build_pedigree(1:3, data = dat)
addtest3 <- additive_genetic_animal(ped3, dat$id)
test_that("The function returns the type expected",{
  expect_output(effect_type(addtest3),"random")
})

## TODO: Test with objects of types: 
##   - genetic
##   - permanent_environmental_competition

