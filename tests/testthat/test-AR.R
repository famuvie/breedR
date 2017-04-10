
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
