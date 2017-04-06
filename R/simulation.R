#' @rdname simulation
#' @name simulation
#' @title Simulation of phenotypes and model components
#' @export breedR.sample.AR breedR.sample.splines breedR.sample.BV
#' @export breedR.sample.phenotype
#' @description These functions allow to draw samples from several models 
#'   (spatial, genetic, competition, etc.) and to combine them to produce a 
#'   simulated phenotype. The resulting dataset can then be fitted with breedR 
#'   to compare the estimations with the true underlying parameters. 
#'   \code{breedR.sample.phenotype} is the main function in the group, as it 
#'   makes use of the rest to simulate a phenotype's components.
#' @param fixed a numeric vector of regression coefficients.
#' @param random a list of random effects specifications, where each element is 
#'   itself a list with elements \code{nlevels} and \code{sigma2}.
#' @param genetic a list with the additive genetic effect specifications. See 
#'   Details.
#' @param spatial a list with the spatial effect specifications. See Details.
#' @param residual.variance is a positive number giving the value of the 
#'   residual variance.
#' @param N number of simulated individuals. If \code{spatial} is specified, 
#'   \code{N} is overrided by the product of \code{spatial$grid.size}. Otherwise
#'   it is required. If \code{genetic} is specified, \code{N} is the size of the
#'   offspring only.
#'   
#' @details The design matrix for the \code{fixed} effects (if given) is a 
#'   column of ones and a matrix of random uniform values in \code{(0, 1)}. 
#'   Therfore, the first element in \code{fixed} gives the overall intercept.
#'   
#'   \code{genetic} is a list with the followng elements:
#'   
#'   \itemize{
#'   
#'   \item \code{model} a character string, either 'add_animal' or 
#'   'competition'. In the former, a single breeding value per individual will 
#'   be simulated, while in the latter \emph{direct} and \emph{competition} 
#'   values are simulated.
#'   
#'   \item \code{Nparents} passed to \code{breedR.sample.pedigree}.
#'   
#'   \item \code{sigma2_a} numeric. For the \code{add_animal} model, the 
#'   variance of the additive genetic effect. For the \code{competition} model, 
#'   the \eqn{2\times 2} covariance matrix of direct and competition genetic 
#'   effects. Passed to \code{breedR.sample.BV} as \code{Sigma}.
#'   
#'   \item \code{check.factorial} passed to \code{breedR.sample.pedigree}
#'   
#'   \item \code{pec} numeric. If present, and only under the \code{competition}
#'   model, it simulates a \emph{Permanent Environmental Competition} effect 
#'   with the given variance.
#'   
#'   \item \code{relations} character. If present and equals \code{half-sibs} it
#'   will generate a pedigree with unknown sires, so that relationships in the 
#'   offsprings are either unrelated or half-sibs are possible. Otherwise, both 
#'   parents are known and full-sibs are also possible.
#'   
#'   }
#'   
#'   Note that only one generation is simulated.
#'   
#'   \code{spatial} is a list with the following elements:
#'   
#'   \itemize{
#'   
#'   \item \code{model} a character string, either 'AR' or 'splines'.
#'   
#'   \item \code{grid.size} a numeric vector of length two with the number of 
#'   rows and columns of trees. Note that the spacing between trees is equal in 
#'   both dimensions.
#'   
#'   \item \code{rho/n.knots} passed to \code{breedR.sample.AR} or to 
#'   \code{breedR.sample.splines} as \code{nkn}.
#'   
#'   \item \code{sigma2_s}  passed to \code{breedR.sample.AR} or to 
#'   \code{breedR.sample.splines} as \code{sigma2}.
#'   
#'   }
#'   
#' @examples
#' 
#' breedR.sample.phenotype(fixed   = c(mu = 10, x = 2),
#'                         random = list(u = list(nlevels = 3,
#'                                                sigma2  = 1)),
#'                         genetic = list(model    = 'add_animal',
#'                                        Nparents = c(10, 10),
#'                                        sigma2_a = 2,
#'                                        check.factorial = FALSE),
#'                         spatial = list(model     = 'AR',
#'                                        grid.size = c(5, 5),
#'                                        rho       = c(.2, .8),
#'                                        sigma2_s  = 1),
#'                         residual.variance = 1)
breedR.sample.phenotype <- function(fixed = NULL,
                                    random = NULL,
                                    genetic = NULL,
                                    spatial = NULL,
                                    residual.variance = 1,
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
    X <- cbind(1,
               matrix(runif(Nfull*(length(fixed) - 1)),
                      nrow = Nfull))
    phenotype <- phenotype + X %*% fixed

    components$X <- as.data.frame(X)
    if( is.null(names(fixed)) ) {
      if( length(fixed) == 1 ) names(fixed) <- 'X'
      else names(fixed) <- paste('X', seq_along(fixed)-1, sep = '')
    }
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
    
    coord <- expand.grid(lapply(spatial$grid.size, seq))
    
    # Randomly distribute the observed trees over the grid
    ord <- sample(Nobs)
    if( !is.null(genetic) ){
      # Make room for founders
      # The first sum(genetic$Nparents) will go out of range and result
      # in coordinates NA
      arrange <- c(Nobs + 1:sum(genetic$Nparents), ord)
    } else arrange <- ord
    
    components <- cbind(components, coord[arrange, ])
    
    if( spatial$model == 'AR') {
      components$spatial  <- breedR.sample.AR(spatial$grid.size,
                                              spatial$rho,
                                              spatial$sigma2_s)[arrange, 'V1']
    } else if( spatial$model == 'splines') {
      if( !exists('n.knots', spatial) ) spatial$n.knots <- NULL
      # this one comes in the right order, because it is already 
      # multiplied by the incidence matrix
      components$spatial  <- rbind(data.frame(V1 =rep(NA, Nfull - Nobs)), 
                               breedR.sample.splines(coord[ord, ],
                                                   spatial$n.knots,
                                                   spatial$sigma2_s))[, 'V1']
    } else stop('Please specify spatial model.')
    
    phenotype <- rowSums(cbind(phenotype, components$spatial), na.rm = TRUE)
  } 
  
  # Genetic
  if( !is.null(genetic) ) {
    
    # Simulate components
    if( is.null(genetic$check.factorial) ) cf <- TRUE
    else cf <- genetic$check.factorial
    
    ped <- suppressWarnings(breedR.sample.pedigree(Nobs,
                                                   genetic$Nparents,
                                                   check.factorial = cf))
    
    # Remove founders without offspring
    if( exists('map', attributes(ped)) ) {
      if ( !cf ) {
        rm.idx <- which(is.na(attr(ped, 'map')))
        components <- components[-rm.idx, ]
        phenotype  <- phenotype[-rm.idx]
        Nfull <- nrow(components)
      } else {
        stop('This should not happen')
      }
    }
    
    # Account for OP pedigrees
    if( exists('relations', genetic) ) {
      if( genetic$relations == 'half-sibs' ) {
        ped@sire <- as.integer(rep(NA, Nfull))
      }
    }
    
    components <- cbind(components, as.data.frame(ped))
    
    ## Include Breeding Values (direct additive and potentially others like comp.)
    ## if sigma2_a is a number, returns a BV column
    ## if it has dimension 2x2, returns a matrix with columnsBV1 and BV2
    components  <- cbind(components, breedR.sample.BV(ped, genetic$sigma2_a))
    
    # Compute the effect of competitors on phenotypes
    if( genetic$model == 'competition' ){
      
      ## Check genetic model
      genetic$pedigree <- ped
      genetic$id <-sum(genetic$Nparents) + 1:Nobs  # index of individuals
      genetic$coord <- coord[ord, ]
      genetic$var.ini <- genetic$sigma2_a
      genetic$response <- 1   # fake response to establish n. traits
      genetic <- do.call(check_genetic, genetic)
      
      ## Build components
      gen_comp <- additive_genetic_competition(
        pedigree    = genetic$pedigree,
        coordinates = genetic$coordinates,
        id          = genetic$id,
        decay       = genetic$competition_decay,
        autofill    = genetic$autofill)
      
      # Incidence matrix (in condensed 8-col format, with neighbour indices)
      Bmat <- renderpf90.matrix( model.matrix(gen_comp))
      Bmat[Bmat==0] <- NA
      
      # Genetic competition values of neighbours
      Cmat <- matrix(components$BV2[Bmat[, 8+1:8]], nrow = Nobs)
      
      # Weighted Neighbour Competition
      components$wnc <- c(rep(NA, Nfull-Nobs),
                          rowSums(Bmat[, 1:8] * Cmat, na.rm = TRUE))
      
      # Permanent Environmental Competition effect
      if (exists('pec', genetic)) {
        components$pec <- c(rep(NA, Nfull-Nobs),
                            rnorm(Nobs, sd = sqrt(genetic$pec$var.ini)))
        Pmat <- matrix(components$pec[Bmat[, 8+1:8]], nrow = Nobs)
        components$wnp <- c(rep(NA, Nfull-Nobs),
                            rowSums(Bmat[, 1:8] * Pmat, na.rm = TRUE))
      }
      
    }
    # Include in the phenotype the corresponding components
    phenotype <- phenotype +
      apply(as.matrix(components[,names(components) %in% c('BV', 'BV1', 'wnc', 'wnp')]),
            1, sum)

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
#'   columns in the field trial
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
  # Note: MASS::mvrnorm only draws samples from the covariance Matrix
  if( !requireNamespace("spam", quietly = TRUE) ) {
    stop("Package spam needed for this simulating AR model. Please install.",
         call. = FALSE)
  }
  ans <- t(spam::rmvnorm.prec(N, Q = Uinv))
  
  return(as.data.frame(ans))
}

#' @rdname simulation
#' @param coord numeric. A two-column matrix(-like) with spatial coordinates.
#' @param nkn numeric. A vector of length two with the number of (inner) knots
#'   in each dimension
#' @details \code{breedR.sample.splines} simulates a two-dimensional spatial 
#'   process as the kronecker product of B-splines processes in each dimension.
breedR.sample.splines <- function(coord, nkn, sigma2, N = 1){
  
  splines_struct <- breedr_splines(coord, nkn)

  Umat <- get_structure(splines_struct)
  
  # Simulated samples of effects
  if( !requireNamespace("MASS", quietly = TRUE) ) {
    stop("Package MASS needed for this simulating splines. Please install.",
         call. = FALSE)
  }
  eff <- t(matrix(MASS::mvrnorm(N, mu = rep(0, dim(Umat)[1]), Sigma = Umat*sigma2),
                  nrow = N))
  
  # Multiply by incidence matrix to get observations at individual level
  ans <- model.matrix(splines_struct) %*% eff
  
  # cast to matrix (rather than possibly Matrix)
  # and then to data.frame (no direct way)
  return(as.data.frame(as.matrix(ans)))
}

#' @rdname simulation
#' @param ped a pedigree object
#' @param Sigma numeric. The additive genetic variance. Either a variance for a 
#'   single additive genetic effect, or a positive-definite matrix with the
#'   covariance structure for a set of correlated genetic effects
#' @details \code{breedR.sample.BV} simulates a set of breeding values (BV)
#'   given a pedigree
breedR.sample.BV <- function(ped, Sigma, N = 1) {
  
  # Precision matrices are more sparse
  Ainv <- pedigreemm::getAInv(ped)
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
#' @param Nparents numeric. Vector of length two. Number of dams and sires to 
#'   randomly mate.
#' @param check.factorial logical. If TRUE (default), it checks whether all the 
#'   possible matings had taken place at least once. If not, it stops with an 
#'   error.
#' @details \code{breedR.sample.pedigree} simulates a one-generation pedigree 
#'   from random mating of independent founders. Note that if 
#'   \code{check.factorial} is \code{FALSE}, you can have some founders removed
#'   from the pedigree.
breedR.sample.pedigree <- function(Nobs, Nparents, check.factorial = TRUE) {
  stopifnot(length(Nparents) == 2)
  if( is.null(names(Nparents)) ) names(Nparents) <- c('mum', 'dad')
  ped.obs <- data.frame(id = 1:Nobs + sum(Nparents),
                        dad = sample(Nparents['dad'],
                                     size = Nobs,
                                     replace = TRUE),
                        mum = sample(Nparents['mum'],
                                     size = Nobs,
                                     replace = TRUE) + Nparents['dad'])
  # Case half-sibs balanced
  if( Nparents['dad'] == 0 & Nobs %% Nparents['mum'] == 0) {
    ped.obs <- data.frame(id = 1:Nobs + sum(Nparents),
                          dad = 0,
                          mum = sample(rep(1:Nparents['mum'], length = Nobs)))
  }
  fullped <- build_pedigree(1:3, data = ped.obs)

  # checks
  if( check.factorial )
    stopifnot(all(xtabs( ~ dad + mum, data = ped.obs) > 0))
  stopifnot(all(check_pedigree(fullped)))
  
  return(fullped)
}


#' @rdname simulation
#' @param dim numeric. Dimension of the effect (e.g. n. of traits)
#' @param var numeric matrix. (Co)variance matrix
#' @param Nlevels numeric. Number of individuals values to sample
#' @param labels character vector of labels for each level.
#' @param N numeric. Number of observations to sample
#' @param vname string. A name for the resulting variables
#' @details \code{breedR.sample.ranef} simulates a random effect with a given
#'   variance.
breedR.sample.ranef <- 
  function(dim, var, Nlevels, labels = NULL, N = Nlevels, vname = 'X') {
  if(!is.null(dim(var))) stopifnot(identical(dim(var), rep(as.integer(dim), 2)))
  if(!is.null(labels)) stopifnot(all.equal(length(labels), Nlevels))
  
  ## Simulate Nlevels correlated vectors of dimension dim
  U <- chol(var)
  values <- matrix(rnorm(dim*Nlevels), ncol = dim) %*% U
  
  ## Sample N observations of the 'factor'
  if (N == Nlevels) {
    idx <- seq_len(Nlevels)
  } else {
    ## either N > Nlevels or N < Nlevels
    idx <- sample(seq_len(Nlevels), N, replace = TRUE)
  }
  ans <- data.frame(values[idx, ])
  
  ## variable names from variance matrix colnames
  if (!is.null(colnames(var))) 
    trait_names <- colnames(var)
  else trait_names <- seq_len(dim)
  varnames <- paste(vname, trait_names, sep = "_")
  
  if (is.null(labels)) {
    ## variable names from variance matrix colnames
    if (!is.null(colnames(var))) 
      labels <- colnames(var)
    else labels <- seq_len(dim)
    names(ans) <- varnames
  } else {
    ## factor labels as an additional variable
    ans <- cbind(factor(labels[idx]), ans)
    names(ans) <- c(vname, varnames)
  }
  
  
  return(ans)
}

