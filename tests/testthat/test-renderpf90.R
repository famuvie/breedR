old.op <- options(warn = -1)  # suppressWarnings
on.exit(options(old.op))

context("renderpf90")

test.inc8x4 <- as(rep(1:4, 2), 'indMatrix')   # An 8x4 incidence matrix
test.cov4x4 <- with(L <- Matrix::tril(matrix(sample(16),4)),
                    Matrix::t(L)%*%L)         # a SPD 4x4 covariance matrix

##Test with an object of type generic
gtest <- generic(incidence  = test.inc8x4, covariance = test.cov4x4)
test_that("The function returns the same number",{
  expect_equal(renderpf90(gtest)$levels,nrow(gtest$structure.matrix))
})

##Test with an object of type splines
x.loc <- 1:100
y.loc <- seq(1000, by = 5, length = 51)
coord <- expand.grid(x.loc, y.loc)
result <- splines(coord)
test_that("The function behaves as expected",{
  expect_equal(tail(renderpf90(result)$levels, 1),nrow(result$structure.matrix))
  expect_true(all(head(renderpf90(result)$levels, -1) ==0))
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
test_that("The function behaves as expected",{
  expect_equal(nrow(renderpf90(comptest)$structure.matrix),nrow(comptest$structure.matrix))
})

##Test with an object of type additive_genetic_competition
ped <- build_pedigree(1:3, data = dat)
addtest <- additive_genetic_competition(ped, coord = coordinates, dat$id, 1, autofill=TRUE)
test_that("The function behaves as expected",{
  expect_equal(nrow(renderpf90(addtest)$structure.matrix),nrow(addtest$structure.matrix))
})

##Test with an object of type additive_genetic
ped2 <- pedigreemm::pedigree(sire = c(NA,NA,1, 1,4,5),
                             dam  = c(NA,NA,2,NA,3,2),
                             label= 1:6)
inc <- cbind(0, 1, diag(4))
addtest2 <- additive_genetic(ped2, inc)
test_that("The function  behaves as expected",{
  expect_equal(nrow(renderpf90(addtest2)$structure.matrix),nrow(addtest2$structure.matrix))
})

##Test with an object of type additive_genetic_animal
dat <- data.frame(id = 1:4,
                  sire = c(11, 11, 2, 3),
                  dam  = c(12, NA, 1, 12))
ped3 <- build_pedigree(1:3, data = dat)
addtest3 <- additive_genetic_animal(ped3, dat$id)
test_that("The function behaves as expected",{
  expect_equal(nrow(renderpf90(addtest3)$structure.matrix),nrow(addtest3$structure.matrix))
})

##Test with an object of type genetic
genetictest <- genetic(ped2, inc,  cov = diag(6))
test_that("The function returns the type expected",{
  expect_equal(nrow(renderpf90(genetictest)$structure.matrix),nrow(genetictest$structure.matrix))
})

##Test with an object of type permanent_environmental_competition
dat_pec <- data.frame(id   = 1:5,
                      sire = c(11, 11, 2, 3, 2),
                      dam  = c(12, NA, 1, 12, 1),
                      x    = c(rep(1:2, times = 2), 3),
                      y    = c(rep(1:2, each = 2), 3))
pectest <- permanent_environmental_competition(coord = dat_pec[, c('x', 'y')], decay = 2)
test_that("The function returns the type expected",{
  expect_equal(nrow(renderpf90(pectest)$structure.matrix),nrow(pectest$structure.matrix))
})

