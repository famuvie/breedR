#' @rdname simulation
#' @name simulation
#' @title Simulation of phenotypes and model components
#' @export breedR.sample.AR breedR.sample.splines breedR.sample.BV 
#'   breedR.sample.phenotype
#' @description These functions allow to draw samples from several models 
#'   (spatial, genetic, competition, etc.) and to combine them to produce a 
#'   simulated phenotype. The resulting dataset can then be fitted with breedR 
#'   to compare the estimations with the true underlying parameters.
#' @param fixed a numeric vector of regression coefficients.
#' @param N numeric. Number of samples to be drawn.
breedR.sample.phenotype <- function(fixed = NULL,
                                    random = NULL,
                                    genetic = NULL,
                                    spatial = NULL,
                                    residual.variance,
                                    N = NULL) {
  
  # Number of observations (measurements)
  if( !is.null(spatial) ) {
    Nobs <- prod(spatial$grid.size)
  } else {
    if( is.null(N) ) stop('Please specify sample size (argument N)\n')
    Nobs <- N
  }
  
  # Number of individuals in the pedigree
  if( !is.null(genetic) ) {
    Nfull <- Nobs + sum(genetic$Nparents)
  } else {
    Nfull <- Nobs
  }

  components <- list()
  phenotype  <- rep(0, Nfull)

  # Fixed
  if( !is.null(fixed) ) {
    X <- matrix(runif(Nfull*length(fixed)),
                nrow = Nfull)
    phenotype <- phenotype + X %*% fixed

    components$X <- as.data.frame(X)
    if( length(fixed) == 1 ) names(fixed) <- 'X'
    else names(fixed) <- paste('X', 1:length(fixed), sep = '')
    names(components$X) <- names(fixed)
  }

  # Random
  if( !is.null(random) ) {
    make.random.single <- function(x, N) {
      lev <- sample(x$nlevels, N, replace = TRUE)
      val <- rnorm(x$nlevels, sd = sqrt(x$sigma2))
      return(factor(val[lev], levels = val))
    }
    
    components$rnd <- as.data.frame(lapply(random,
                                           make.random.single,
                                           N = Nfull))
    phenotype <- phenotype + apply(sapply(components$rnd,
                                          function(x) as.numeric(levels(x))[x]),
                                   1, sum)
  }
  
  # Spatial
  if( !is.null(spatial) ) {
    
    coord <- expand.grid(sapply(spatial$grid.size, seq))
    
    # Randomly distribute the observed trees over the grid
    ord <- sample(Nobs)
    if( !is.null(genetic) ){
      # Make room for founders
      arrange <- c(Nobs + 1:sum(genetic$Nparents), ord)
    } else arrange <- ord
    
    components <- cbind(components, coord[arrange, ])
    components$AR  <- breedR.sample.AR(spatial$grid.size,
                                       spatial$rho,
                                       spatial$sigma2_s)[arrange, ]
    
    phenotype <- rowSums(cbind(phenotype, components$AR), na.rm = TRUE)
  } 
  
  # Genetic
  if( !is.null(genetic) ) {
    
    # Simulate components
    if( is.null(genetic$check.factorial) ) cf <- TRUE
    else cf <- genetic$check.factorial
    
    ped <- breedR.sample.pedigree(Nobs, genetic$Nparents,
                                  check.factorial = cf)
    components <- cbind(components, as.data.frame(ped))
    components$BV  <- breedR.sample.BV(ped, genetic$sigma2_a)
    
    phenotype <- phenotype + as.numeric(components$BV)
  } else genetic$Nparents = 0

  # Residual
  components$resid <- rnorm(Nfull, sd = sqrt(residual.variance))
  phenotype <- phenotype + components$resid
  
  # Phenotype
  components$phenotype <- phenotype
  
  # Preferred order of the columns for the first elements
  pref.ord <- c('self', 'sire', 'dam', 'row', 'col')
  reord <- na.omit(match(pref.ord, names(components)))
  if( length(reord) > 0 ) {
    components <- c(components[, reord], components[, -reord])
  }
  
  dat <- do.call('data.frame', components)

  return(dat)
}

#' @rdname simulation
#' @param size numeric. A vector of length two with the number of rows and
#'   columns
#' @param rho numeric. A vector of length two with the autocorrelation 
#'   parameters for the row and column autoregressive processes
#' @param sigma2 numeric. The marginal variance
#' @details \code{breedR.sample.AR} simulates a two-dimensional spatial process 
#'   as the kronecker product of first-order autoregressive processes in each 
#'   dimension.
breedR.sample.AR <- function(size, rho, sigma2, N = 1){

  # Precision matrices of 1D-AR processes (unscaled)
  Q1d <- mapply(build.AR1d, size, rho)
  
  # Precision matrix for the AR1(rho_x) x AR1(rho_y) process
  # when locations are stacked following the standards of R:
  # by columns: first vary x and then y
  # (1, 1), (2, 1), ..., (n_x, 1), (1, 2), ..., (n_x, 2), ...
  # Only the lower triangle
  Uinv <- kronecker(Q1d[[2]], Q1d[[1]])
  
  scaling <- prod(1/(1-rho**2))/sigma2
  Uinv <- Uinv * scaling
  
  # Simulated samples 
  # Note: MASS:mvrnorm only draws samples from the covariance Matrix
  ans <- t(spam::rmvnorm.prec(N, Q = Uinv))
  
  if(N == 1) colnames(ans) <- 'AR'
  else colnames(ans) <- paste('AR', 1:N, sep = '')
  
  return(as.data.frame(ans))
}

#' @rdname simulation
#' @param size numeric. A vector of length two with the number of rows and 
#'   columns
#' @param nkn numeric. A vector of length two with the number of knots in each 
#'   dimension parameters for the row and column autoregressive processes
#' @param sigma2 numeric. The marginal variance
#' @details \code{breedR.sample.splines} simulates a two-dimensional spatial
#'   process as the kronecker product of B-splines processes in each dimension.
breedR.sample.splines <- function(size, nkn, sigma2, N = 1){
  
  
  return(NULL)
}

#' @rdname simulation
#' @param ped a pedigree object
#' @param Sigma numeric. The additive genetic variance. Either a variance for a
#'   single additive genetic effect, or a matrix with the covariance structure for a set of
#'   correlated genetic effects
#' @details \code{breedR.sample.BV} simulates a set of breeding values (BV) given a pedigree
breedR.sample.BV <- function(ped, Sigma, N = 1) {
  
  # Precision matrices are more sparse
  Ainv <- getAInv(ped)
  Q <- kronecker(solve(Sigma), Ainv)
  
  # number of genetic effects
  neff <- dim(as.matrix(Sigma))[1]
  
  # A dim(Ainv)[1] x dim(Sigma)[1] vector with simulated BV
  ans <- matrix(spam::rmvnorm.prec(N, Q = Q),
                nrow = dim(Ainv)[1])
  if(neff == 1) colnames(ans) <- 'BV'
  else colnames(ans) <- paste('BV', 1:neff, sep ='')
  
  return(ans)
}


#' @rdname simulation
#' @param Nobs numeric. Number of individuals to sample
#' @param Nparents numeric. Vector of length two. Number of fathers and mothers
#'   to randomly mate.
#' @param check.factorial. logical. If TRUE, checks whether all the possible
#'   matings had taken place at least once.
#' @details \code{breedR.sample.pedigree} simulates a one-generation pedigree
#'   from random mating of independent founders
breedR.sample.pedigree <- function(Nobs, Nparents, check.factorial = TRUE) {
  stopifnot(length(Nparents) == 2)
  names(Nparents) <- c('dad', 'mum')
  ped.obs <- data.frame(id = 1:Nobs + sum(Nparents),
                        dad = sample(Nparents['dad'],
                                     size = Nobs,
                                     replace = TRUE),
                        mum = sample(Nparents['mum'],
                                     size = Nobs,
                                     replace = TRUE) + Nparents['dad'])
  fullped <- build_pedigree(1:3, data = ped.obs)
  
  # checks
  if( check.factorial )
    stopifnot(all(xtabs( ~ dad + mum, data = ped.obs) > 0))
  stopifnot(all(check_pedigree(fullped)))
  
  return(fullped)
}
