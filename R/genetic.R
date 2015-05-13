#' Build a genetic model
#' 
#' Check conformity of arguments and return a \code{genetic} object.
#' 
#' This is a virtual class. No objects are expected to be created directly.
#' 
#' @param pedigree object of class 'pedigree'
#' @inheritParams random
#' @return A list with elements \code{pedigree}, \code{incidence.matrix},
#'   \code{structure.matrix} and \code{structure.type}, which is a string
#'   indicating either \code{covariance} or \code{precision}.
genetic <- function(pedigree, incidence, covariance, precision) {
  
  ## checks
  stopifnot(!missing(pedigree))
  if (nrow(as.data.frame(pedigree)) != ncol(incidence))
    stop('The incidence matrix should have as many columns as individuals in the pedigree.')
  
  ## Build the random effect, and further specify the genetic class
  random.call <- mc <- match.call()
  arg.list <- as.list(mc)[-1]
  
  random.call[[1]] <- as.symbol('random')
  ans <- eval(random.call[-match('pedigree', names(random.call))],
              parent.frame())

  ans$pedigree <- pedigree
  class(ans) <- c('genetic', class(ans))
  
  return(ans)
}


#' Build an additive_genetic model
#' 
#' Check conformity of arguments and return a \code{additive_genetic} object.
#' 
#' @param pedigree object of class 'pedigree'
#' @inheritParams random
#' @return A list with elements \code{pedigree}, \code{incidence.matrix},
#'   \code{structure.matrix} and \code{structure.type}, which is a string
#'   indicating either \code{covariance} or \code{precision}.
#' @examples 
#' ped <- pedigreemm::pedigree(sire = c(NA,NA,1, 1,4,5),
#'                             dam  = c(NA,NA,2,NA,3,2),
#'                             label= 1:6)
#' inc <- cbind(0, 0, diag(4))
#' additive_genetic(ped, inc)
additive_genetic <- function(pedigree, incidence) {
  
  
  ## Build the genetic effect, and further specify the additive_genetic class
  relationship.matrix <- pedigreemm::getA(pedigree)
  
  ## NOTE:
  ## We could as well use pedigreemm::getAInv(pedigree)
  ## to construct the precision matrix.
  
  ans <- genetic(pedigree, incidence, covariance = relationship.matrix)

  class(ans) <- c('additive_genetic', class(ans))
  
  return(ans)
}


#' Build an additive-genetic animal model
#' 
#' Given a pedigree, and an index vector of observations, build and 
#' \code{additive_genetic_animal} model.
#' 
#' \code{idx} must hold the index of observed individuals in the original 
#' codification. If recoding took place when building the pedigree, this
#' function will convert the codes internally.
#' 
#' @param idx integer vector of observed individuals (in the original codif.)
#' @inheritParams additive_genetic
#'   
#' @return A list with elements \code{pedigree}, \code{incidence.matrix}, 
#'   \code{structure.matrix} and \code{structure.type}, which is a string 
#'   indicating either \code{covariance} or \code{precision}.
#' @examples 
#' dat <- data.frame(id = 1:4,
#'                   sire = c(11, 11, 2, 3),
#'                   dam  = c(12, NA, 1, 12))
#' ped <- build_pedigree(1:3, data = dat)
#' additive_genetic_animal(ped, dat$id)
additive_genetic_animal <- function(pedigree, idx) {
  
  ## Checks
  stopifnot(is.numeric(idx))
  stopifnot(length(idx) < nrow(as.data.frame(pedigree)))
  
  ## Incidence matrix
  ## It is possible that the pedigree has been recoded/reordered
  ## In that case, we need to recode the data file id codes as well
  if (!is.null(attr(pedigree, 'map')))
    idx <- attr(pedigree, 'map')[idx]
  
  inc.mat <- as(idx, 'indMatrix')
  
  ans <- additive_genetic(pedigree, inc.mat)
  class(ans) <- c('additive_genetic_animal', class(ans))
  return(ans)
}



#' Build an additive-genetic competition model
#' 
#' Return incidence and structure for a \code{additive_genetic_competition}
#' model, given the pedigree, the spatial coordinates and codes of the
#' observations and the competition decay parameter.
#' 
#' \code{id} must hold the codes of observed individuals in the original 
#' codification. If recoding took place when building the pedigree, this 
#' function will handle the codes internally.
#' 
#' @param id integer vector of numeric codes for observed individuals
#' @inheritParams additive_genetic
#' @inheritParams competition
#' @inheritParams build_grid
#'   
#' @return A list with elements \code{pedigree}, \code{incidence.matrix}, 
#'   \code{structure.matrix} and \code{structure.type}, which is a string 
#'   indicating either \code{covariance} or \code{precision}.
#' @examples 
#' dat <- data.frame(id   = 1:5,
#'                   sire = c(11, 11, 2, 3, 2),
#'                   dam  = c(12, NA, 1, 12, 1),
#'                   x    = c(rep(1:2, times = 2), 3),
#'                   y    = c(rep(1:2, each = 2), 3))
#' ped <- build_pedigree(1:3, data = dat)
#' additive_genetic_competition(ped, coord = dat[, c('x', 'y')], dat$id, 2)
additive_genetic_competition <- function(pedigree,
                                         coordinates,
                                         id,
                                         decay,
                                         autofill = TRUE) {
  
  ## Checks
  stopifnot(is.numeric(id))
  stopifnot( (n <- length(id)) < (p <- nrow(as.data.frame(pedigree))))
  stopifnot(length(id) == nrow(coordinates))
  
  ## additive_genetic_competition inherits from additive_genetic
  ## and from competition simultaneously.
  ## As S3 classes do not support multiple inheritance, simulate it
  ## manually by constructing separate competition and additive_genetic
  ## objects with dummy covariance and incidence matrices respectively
  ## and then composing manually the random effect
  cov.dummy <- diag(n)
  comp.aux <- competition(coordinates = coordinates, 
                          covariance  = cov.dummy,
                          decay       = decay,
                          autofill    = autofill)
  
  inc.dummy <- diag(p)
  ag.aux <- additive_genetic(pedigree, inc.dummy)
  
  ## the internal codes of the observed individuals are the indices
  ## of the corresponding levels of the random effect
  id.internal  <- attr(pedigree, 'map')[id]
  
  ## the incidence matrix has to be exapanded with 0
  ## in the columns corresponding to unobserved individuals (e.g. founders)
  ## furthermore, the order of the columns must fit
  ## the internal representation
  inc.mat <- matrix(0, n, p)
  inc.mat[, id.internal] <- as.matrix(comp.aux$incidence.matrix)
  
  random.args <- structure(list(inc.mat, ag.aux$structure.matrix),
                           names = c('incidence', ag.aux$structure.type))
  
  ans <- do.call('random', random.args)
  ans$pedigree <- pedigree
  ans$coordinates <- coordinates
  
  ## Include inherited classes from the additive_genetic object 
  ## in the first place
  class(ans) <- c('additive_genetic_competition',
                  setdiff(class(ag.aux), class(ans)),
                  setdiff(class(comp.aux), class(ans)),
                  class(ans))
  
  return(ans)
}



# Build additive direct and competition models
build.genetic.model <- function(genetic) {
  
  # Checks
  
  
  # Incidence matrix
  # It is possible that the pedigree has been recoded/reordered
  # In that case, we need to recode the data file id codes as well
  if( !is.null(attr(genetic$pedigree, 'map')) )
    idx <- attr(genetic$pedigree, 'map')[genetic$id]
  else
    idx <- genetic$id
  
  B = as.matrix(idx, ncol = 1)
  
  ans <- list()
  
  # In a competition model, there are further columns to be added to the right
  if( genetic$model == 'competition' ) {
    
    # Original coordinates
    coord0 <- as.data.frame(sapply(genetic$coord, as.numeric))
    
    # lattice of spatial locations
    # possibly with automatic filling of empty rows or columns
    pos <- loc_grid(coord0, genetic$autofill)
    
    # The coordinates as factors allows to find easily the number of 
    # different x and y positions, and the ordering
    coord <- as.data.frame(mapply(factor,
                                  coord0,
                                  pos,
                                  SIMPLIFY = FALSE))
    
    # Number of different locations in rows and cols
    pos.length <- sapply(pos, length)
    
    # Spacing between trees in rows and columns
    pos.step <- sapply(pos, function(x) diff(x)[1])
    
    # Map data coordinates with corresponding index of the Q matrix
    matrix2vec <- function(x, nx = pos.length[1], ny = pos.length[2]) {
      map <- matrix(1:(nx*ny), nx, ny)
      return(apply(x, 1, function(y) map[y[1], y[2]]))
    }
    ord <- matrix2vec(sapply(coord, as.integer))
    
    # TODO: should I abstract this grid and coordinates procedures?
    
    # Check: is this a regular grid?
    # if regular, n_x \times n_y ~ n_obs
    # if irregular, n_x \times n_y ~ n_obs^2
    # We assume it is regular if
    # n_x \times n_y < n_obs + (n_obs^2 - n_obs)/4
    if( prod(pos.length) > nrow(coord)*(nrow(coord) + 3)/4 )
      stop('The competition model can only be fitted to regular grids.')
    
    # Check: the grid should be actually 2d
    if( !all(pos.length > 2) )
      stop('Are you kidding? This is a line!')
    
    # Each individual has (up to) eight neighbours.
    # Moving clockwise from North, this matrix stores the id
    # of the corresponding neighbourg, or zeroes if missing
    dirs <- c('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW')
    
    # Distances in each direction
    diag.dist <- sqrt(sum(pos.step**2))
    dir.dist <- rep(c(pos.step[1], diag.dist, pos.step[2], diag.dist), 2)
    
    # Distance weighting in each direction
    wdir.dist <- 1/dir.dist**genetic$competition_decay
    
    # Matrix of trees idx in their respective position
    idx.mat <- matrix(NA, pos.length[1], pos.length[2])
    idx.mat[ord] <- idx
    
    # idx of neighbour indices at each direction
    Bneigh <- neighbours.at(idx.mat, dirs)[ord, ]
    
    # IC of neighbouring trees at each direction
    # normalize to make all coefficients-squared add up to one
    BIC <- sapply(seq_along(dirs),
                  function(dir) ifelse(is.na(Bneigh[, dir]), NA, wdir.dist[dir]))
    BIC <- BIC / sqrt(apply(BIC**2, 1, sum, na.rm = TRUE))
    
    # check
    stopifnot(all(sapply(apply(BIC**2, 1, sum, na.rm = TRUE), all.equal, 1)))
    
    
      
    # cbind the Intensity of Competition and the neighbours idx
    B <- cbind(B, BIC, Bneigh)

    ans <- list(param = genetic$competition_decay,
                coord = coord0)
  }

  B[is.na(B)] <- 0
  return(c(ans, list(B = B)))
}
