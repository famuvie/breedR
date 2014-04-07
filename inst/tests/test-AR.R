#### Build small testbeds ####
build.testbed <- function(corner = c(0, 0), size, treesep = c(1, 1), beta){
  n = size[1] * size[2]
  s.mat = matrix(NA, nrow = size[1], ncol = size[2])
  j = 1:size[2]
  for(i in 1:size[1])
    s.mat[i,j] = 0.1*(i+2*j)
  ## a covariate
  set.seed(2)
  z.mat = matrix(runif(n), size[1], size[2])
  ## noise
  set.seed(2)
  noise.mat = matrix(rnorm(n, sd = 0.3), size[1], size[2])
  ## make simulated data
  y.mat = beta * z.mat + s.mat + noise.mat
  ## build final dataset  
  dat <- data.frame(i = rep(seq(corner[1], by = treesep[1], length = size[1]),
                            times = size[2]),
                    j = rep(seq(corner[2], by = treesep[2], length = size[2]),
                            each = size[1]),
                    z = as.vector(z.mat),
                    y = as.vector(y.mat),
                    true.s = as.vector(s.mat))
  return(dat)
}

beta = 0.5
datlist <- list(# small square regular grid
  small.sq.reg = build.testbed(corner = c(0, 0),
                               size = c(5, 5),
                               treesep = c(1, 1),
                               beta = beta),
  # small rectangular grid with different spacings and coordinates
  small.rect.irr = build.testbed(corner = c(134, 77),
                                 size = c(5, 7),
                                 treesep = c(3, 4),
                                 beta = beta))

# triangular configuration
datlist <- c(datlist,
             triang = list(datlist[[1]][which(as.vector(tril(matrix(TRUE, 5, 5)))),]))



#### Context: AR models with diffferent arrangements of trees ####
context("AR models with diffferent arrangements of trees")


# Fit models both with EM and AI-REML
run.model <- function(dat, method) {
  res = try(suppressWarnings(remlf90(fixed = y ~ 1 + z, 
                                     spatial = list(model = 'AR',
                                                    coord = dat[, 1:2],
                                                    rho = c(.9, .9)),
                                     data = dat,
                                     method = method)),
                             silent = TRUE)
  return(list(dat = dat,
              method = method,
              res = res))
}

reslist <- c(llply(datlist, run.model, method = 'em'),
             llply(datlist, run.model, method = 'ai'))
             
# Check results
# summary(reslist[[1]])
# res <- reslist[[1]]
# dat <- datlist[[1]]
# require(plyr)
check.result <- function(m, datlabel, debug.plot = FALSE) {
  test_that(paste("AR model runs OK with dataset", datlabel, "and method", m$method), {
    expect_true(!inherits(m$res, 'try-error'))
  })
  
  if( !inherits(m$res, 'try-error') ){
    # Mean Square Error for the spatial effect
    if(debug.plot) {
      print(qplot(as.vector(m$dat$true.s), fixef(m$res)$Intercept$value + m$res$spatial$fit$z) +
              geom_abline(int = 0, sl = 1))
    }
    mse <- mean((as.vector(m$dat$true.s) - fixef(m$res)$Intercept$value + m$res$spatial$fit$z)^2)
    test_that(paste("MSE of the spatial effect estimation is reasonable for dataset",
                    datlabel, "and method", m$method), {
      expect_that(mse, is_less_than(1))
    })
    
    # Estimate of the linear coefficient
    beta.e <- beta - fixef(m$res)$z$value
    test_that(paste("The linear coefficient is estimated within 3 se for dataset",
              datlabel, "and method", m$method), {
      expect_that(abs(beta.e), is_less_than(3*fixef(m$res)$z$s.e.))
    })
  }
}

for(i in 1:length(reslist)) 
  check.result(reslist[[i]], names(reslist)[i], debug.plot = FALSE)



#### Context: selection of autoregressive parameters ####
context("Selection of autoregressive parameters")

res.unset <- try(remlf90(fixed = y ~ z, 
                         spatial = list(model = 'AR',
                                        coord = datlist[[1]][, 1:2]),
                         data = datlist[[1]]),
                 silent = TRUE)

test_that("if rho unset, remlf90 tries a grid of combinations", {
  # remlf90() returns an evaluation grid
  expect_that(exists('rho', as.environment(res.unset)), is_true())
  # the evaluation grid returns the loglikelihood for each default combination
  expect_that(all(complete.cases(res.unset$rho$loglik)), is_true())
})


gridlist <- list(expand.grid(seq(80, 90, 5), c(87, 93))/100,
                 expand.grid(seq(80, 90, 5), NA)/100,
                 expand.grid(NA, c(87, 93))/100)
reslist.spec <- lapply(gridlist, function(g)
  try(remlf90(fixed = y ~ z, 
              spatial = list(model = 'AR',
                             coord = datlist[[1]][, 1:2],
                             rho = g),
              data = datlist[[1]]),
      silent = TRUE))

test_that("the user can specify a full or partial grid of combinations", {
  
  for(i in 1:length(reslist.spec)) {
    res <- reslist.spec[[i]]
    grid <- gridlist[[i]]
    
    # remlf90() returns an evaluation grid
    expect_that(exists('rho', as.environment(res)), is_true())
    
    # The evaluation grid conforms to the user specification
    get_levels <- function(levels) {
      if(all(is.na(levels))) return(breedR.getOption('ar.eval'))
      else return(levels)
    }
    eval.grid <- expand.grid(lapply(lapply(grid, unique), get_levels),
                             KEEP.OUT.ATTRS = FALSE)
    names(eval.grid) <- names(res$rho)[1:2]
    expect_identical(eval.grid,
                     res$rho[, 1:2])

    # the evaluation grid returns the loglikelihood for each combination specified
    expect_that(all(complete.cases(res$rho$loglik)),
                is_true())
  }
})


# # Debug
# image(s.mat)
# image(matrix(res.bR$spatial$fit$z, nrow, ncol))
# qplot(as.vector(s.mat), res.bR$spatial$fit$z) + geom_abline(int = 0, sl = 1, col = 'darkgray')
# summary(res.bR)

