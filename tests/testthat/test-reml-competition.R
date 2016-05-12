old.op <- options(warn = -1,  # suppressWarnings
                  show.error.messages = FALSE)  # silent try
on.exit(options(old.op))

### For testing competition, we perform a simulation excercise ###

set.seed(12345)
grid.size <- c(x = 20, y = 25)   # x: columns, y: rows
# dist.rc <- c(x = 3, y = 5)
coord <- expand.grid(sapply(grid.size, seq))
Nobs <- prod(grid.size)
Nparents <- c(mum = 20, dad = 20)
rho <- -.7  # genetic correlation between additive and competitive values
            # covariance = -.7 * sqrt(2) = -0.98
sigma2_a <- 2   # additive genetic variance
sigma2_c <- 1   # competitive genetic variance
sigma2_p <- .5  # permanent environment effect of competition
sigma2   <- .5   # residual variance

ped.obs <- data.frame(id = 1:Nobs + sum(Nparents),
                      mum = sample(Nparents['mum'],
                                   size = Nobs,
                                   replace = TRUE),
                      dad = sample(Nparents['dad'],
                                   size = Nobs,
                                   replace = TRUE) + Nparents['mum'])
fullped <- build_pedigree(1:3, data = ped.obs)

# checks
# stopifnot(all(xtabs( ~ mum + dad, data = ped.obs) > 0))
stopifnot(all(check_pedigree(fullped)))

# Additive matrix
# Precision matrices are more sparse
  # # Workaround with package pedigree
  # makeAinv(as.data.frame(fullped))
  # Ai <- read.table("Ainv.txt")
  # nInd <- nrow(as.data.frame(fullped))
  # Ainv <- matrix(0,nrow = nInd,ncol = nInd)
  # Ainv[as.matrix(Ai[,1:2])] <- Ai[,3]
  # dd <- diag(Ainv)
  # Ainv <- Ainv + t(Ainv)
  # diag(Ainv) <- dd
Ainv <- pedigreemm::getAInv(fullped)

# Covariance matrix for the genetic effects (a, c)
S_ac <- matrix(c(sigma2_a, rho*sqrt(sigma2_a*sigma2_c),
                 rho*sqrt(sigma2_a*sigma2_c), sigma2_c), 2, 2)
Q <- kronecker(solve(S_ac), Ainv)

# Simulate genetic effects
ped <- as.data.frame(fullped)
gen_a_c <- matrix(spam::rmvnorm.prec(n = 1,
                                     Q = Q),
                  ncol = 2,
                  dimnames = list(NULL, c('a', 'c')))

# checks
# plot(gen_a_c, pch = 19)
# cor(gen_a_c)

# Randomly distribute the trees over the grid
dat <- data.frame(coord[sample(Nobs),],
                  ped.obs,
                  tail(gen_a_c, Nobs),
                  pec = rnorm(Nobs, sd = sqrt(sigma2_p)),
                  e = rnorm(Nobs, sd = sqrt(sigma2)))

# check
# BVs must be more similar for members of the same family
# with(dat, qplot(paste(mum, dad, sep='_'), c))

# Simulate phenotype
# Each tree is affected by the competition value of its neighbours

# Compute IFCs and corresponding neighbours. Assume equal row/col spacing and
# competition intensity decreasing with inverse distance


# neighbours and weighted coefficients of IC
X <- local({
  rect.dist = 1
  diag.dist = sqrt(2)

  # It is convenient to work with matrices representing the spatial arrangement
  ord <- order(dat$x, dat$y)
  matlst <- lapply(dat[ord, c('id', 'a', 'c', 'e', 'pec')],
                   function(x) matrix(x, nrow = grid.size['y']))
  
  rect <- breedR:::neighbours.at.list(matlst, c('N', 'S', 'E', 'W'))
  diag <- breedR:::neighbours.at.list(matlst, c('NE', 'SE', 'SW', 'NW'))

  dat <- c(rect$id,
           diag$id,
           ifelse(is.na(rect$id), NA, 1/rect.dist),
           ifelse(is.na(diag$id), NA, 1/diag.dist),
           rect$c,
           diag$c,
           rect$pec,
           diag$pec)

  # Four-dimensional array
  # two first dimensions are spatial
  # third dimension is direction of neighbourhood: N, S, ..., NW (8)
  # fourth dimension is 
  # 1 = neighbour idx
  # 2 = neighbour IFC
  # 3 = neighbour c
  # 4 = neighbour pec
  x <- array(dat, dim = c(rev(grid.size), 8, 4))
  
  # normalize to make all squared-coefficients add up to one throughout dim. 3
  normalizing.constant = apply(x[,,,2]**2, 1:2, sum, na.rm = TRUE)
  x[,,,2] <- x[,,,2] / as.vector(sqrt(normalizing.constant))
  
  # check
  stopifnot(all(sapply(apply(x[,,,2]**2, 1:2, sum, na.rm = TRUE), all.equal, 1)))
  
  # result in tabular form
  # first eight cols are neighbour idx, last eight are coefs that add up to 1
  res <- data.frame(as.vector(matlst$id),
                    array(x, dim = c(prod(grid.size), 32)))
  colnames(res) <- c('idx',
                     paste('n', 1:8, sep = ''),
                     paste('ifc', 1:8, sep = ''),
                     paste('c', 1:8, sep = ''),
                     paste('pec', 1:8, sep = ''))
  rownames(res) <- NULL
  
  # check
  stopifnot(all(sapply(apply(res[,10:17]**2, 1, sum, na.rm = TRUE), all.equal, 1)))
  
  # return results in the original order of the data frame
  res[order(ord),]
  })

# Contribution to phenotype of neighbouring genetic-competition effects
dat$wnc <- rowSums(X[,10:17] * X[,17+1:8], na.rm = TRUE)

dat$pec <- rowSums(X[,10:17] * X[,25+1:8], na.rm = TRUE)

# Simulated phenotype (with or without pec)
dat <- transform(dat, z = a + wnc + pec + e)
# dat <- transform(dat, z = a + wnc + e)



#### Fitting the competition model with remlf90 ################################
context('Fitting competition models')

expect_error(
  res <- suppressMessages(
    remlf90(
      fixed  = z ~ 1,
      genetic = list(model = c('comp'), 
                     pedigree = dat[, c('id', 'mum', 'dad')],
                     id = 'id',
                     coord = dat[, c('x', 'y')],
                     competition_decay = 1,
                     pec = list(present = TRUE)), 
      data = dat,
      method = 'em',
      debug = F)
  ),
  NA)


# ggplot2::qplot(dat$z - dat$e, fitted(res)) + 
#   ggplot2::geom_abline(intercept = 0, slope = 1, col = 'darkgray')

#### Context: Extraction of results from competition model #####################
context("Extraction of results from competition model")


n.fixed   <- 1
nlevels.fixed <- 1
n.bvs <- nrow(as.data.frame(fullped))  # one set for direct, another for comp.
n.pec <- nrow(dat)

test_that("coef() gets a named vector of coefficients", {
  expect_is(coef(res), 'numeric')
  expect_equal(length(coef(res)), nlevels.fixed + 2*n.bvs + n.pec)
  expect_named(coef(res))
})

test_that("ExtractAIC() gets one number", {
  expect_is(extractAIC(res), 'numeric')
  expect_equal(length(extractAIC(res)), 1)
})

test_that("fitted() gets a vector of length N", {
  expect_is(fitted(res), 'numeric')
  expect_equal(length(fitted(res)), Nobs)
})

test_that("fixef() gets a named list of data.frames with estimated values and s.e.", {
  x <- fixef(res)
  expect_is(x, 'list')
  expect_named(x)
  expect_equal(length(x), n.fixed)
  for (f in x) {
    expect_is(f, 'data.frame')
    expect_named(f, c('value', 's.e.'))
  }
})

test_that("get_pedigree() returns the given pedigree", {
  expect_identical(get_pedigree(res), fullped)
})

test_that("logLik() gets an object of class logLik", {
  expect_is(logLik(res), 'logLik')
})

test_that("model.frame() gets an Nx2 data.frame with a 'terms' attribute", {
  x <- model.frame(res)
  expect_is(x, 'data.frame')
  expect_is(terms(x), 'terms')
  expect_equal(dim(x), c(Nobs, n.fixed + 1))
})

test_that("model.matrix() gets a named list of incidence matrices", {
  x <- model.matrix(res)
  expect_is(x, 'list')
  expect_named(x, get_efnames(res$effects))
  expect_equal(dim(x$Intercept), c(Nobs, nlevels.fixed))
  for (m in x$random) {
    expect_is(m, 'sparseMatrix')
  }
  expect_equal(dim(x$genetic_direct), c(Nobs, n.bvs)) 
  expect_equal(dim(x$genetic_competition), c(Nobs, n.bvs)) 
  expect_equal(dim(x$pec), c(Nobs, Nobs)) 
})

test_that("nobs() gets the number of observations", {
  expect_equal(nobs(res), Nobs)
})

test_that("plot(, type = *) returns ggplot objects after providing coords", {
  ## Even when there is no spatial effect,
  ## the coordinates must be recovered from the competition
  expect_is(plot(res, type = 'phenotype'), 'ggplot')
  expect_is(plot(res, type = 'fitted'), 'ggplot')
  expect_is(plot(res, type = 'residuals'), 'ggplot')
  
  ## But still get errors for the absent spatial components
  expect_error(plot(res, type = 'spatial'), 'no spatial effect')
  expect_error(plot(res, type = 'fullspatial'), 'no spatial effect')
})

test_that("print() shows some basic information", {
  ## Not very informative currently...
  expect_output(print(res), 'Data')
})

test_that("ranef() gets a ranef.breedR object with random effect BLUPs and their s.e.", {
  x <- ranef(res)
  expect_is(x, 'ranef.breedR')
  expect_equal(length(x), 3)
  expect_named(x, c('genetic_direct', 'genetic_competition', 'pec'))
  
  for (y in x) {
    expect_is(y, 'numeric')
    #expect_equal(length(y), n.bvs) # pec is of length 500 
    expect_false(is.null(xse <- attr(y, 'se')))
    
    expect_is(xse, 'numeric')
    #expect_equal(length(xse), n.bvs)
  }
})

test_that("residuals() gets a vector of length N", {
  rsd <- residuals(res)
  expect_is(rsd, 'numeric')
  expect_equal(length(rsd), Nobs)
})

test_that("summary() shows summary information", {
  expect_output(print(summary(res)), 'Variance components')
})

test_that("vcov() gets the covariance matrix of the genetic components of the observations", {

  for (ef in paste('genetic', c('direct', 'competition'), sep = '_')) {
  x <- try(vcov(res, effect = 'genetic_direct'))
  expect_false(inherits(x, 'try-error'))
  expect_is(x, 'sparseMatrix')
  expect_equal(dim(x), rep(Nobs, 2))
  }
})

