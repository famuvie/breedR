#' Build a virtual competition model
#' 
#' Given the coordinates of a set of observations, a decay parameter and a 
#' structure matrix, compute the incidence matrix of competition, and return a 
#' random effect with the given structure.
#' 
#' Exactly one of \code{covariance} and \code{precision} must be specified.
#' 
#' The competition model attributes to each individual a random effect of 
#' competition with variance \eqn{\sigma_{a_c}^2}, which impacts the phenotype 
#' of the neighbours rather than its own phenotype.
#' 
#' Conversely, the effect of the competition over one's phenotype is given by 
#' the additive-genetic competition effects of the neighbours, weighted by the 
#' relative distances. If \eqn{\alpha} is the decay parameter and \eqn{a_c} is 
#' the random competition effect of a neighbour at distance \eqn{d}, then the
#' Weighted Neighbour Competition effect over one's phenotype is given by
#' \deqn{wnc = \sum_{\mathrm{neighbours}} k d^{-\alpha} a_c,}{% wnc =
#' \sum_{neighbours} k (1/d)^\alpha a_c,} where \eqn{k} is a normalizing
#' constant which makes \eqn{Var(wnc) = \sigma_{a_c}^2} and independent of the
#' number of neighbours.
#' 
#' @param decay numeric. The positive value of the decay parameter \eqn{\alpha}.
#'   Typically 1 or 2. See Details.
#' @inheritParams build_grid
#' @inheritParams random
#'   
#' @return An object inheriting from \code{spatial}.
competition <- function(coordinates,
                        covariance,
                        precision,
                        decay,
                        autofill = TRUE) {
  
  mc <- match.call()
  
  ## Checks
  stopifnot(ncol(coordinates) == 2)
  
  coordinates <- as.data.frame(coordinates)
  
  ## Incidence matrix
  grid <- build_grid(coordinates, autofill)
  
  # Check: is this a regular grid?
  # if regular, n_x \times n_y ~ n_obs
  # if irregular, n_x \times n_y ~ n_obs^2
  # We assume it is regular if
  # n_x \times n_y < n_obs + (n_obs^2 - n_obs)/4
  if( prod(grid$length) > nrow(coordinates)*(nrow(coordinates) + 3)/4 )
    stop('The competition model can only be fitted to regular grids.')
  
  # Check: the grid should be actually 2d
  if( !all(grid$length > 2) )
    stop('Are you kidding? This is a line!')
  
  # Each individual has (up to) eight neighbours.
  # Moving clockwise from North, this matrix stores the id
  # of the corresponding neighbourg, or zeroes if missing
  dirs <- c('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW')
  
  # Distances in each direction
  diag.dist <- sqrt(sum(grid$step**2))
  dir.dist <- structure(rep(c(grid$step[1], diag.dist, grid$step[2], diag.dist),
                            2),
                        names = dirs)
  
  # Distance weighting in each direction
  wdir.dist <- 1/dir.dist**decay
  
  # Matrix of trees idx in their respective position
  idx.mat <- matrix(NA, grid$length[1], grid$length[2])
  idx.mat[grid$idx] <- 1:nrow(coordinates)
  
  # idx of neighbour indices at each direction
  Bneigh <- neighbours.at(idx.mat, dirs)[grid$idx, ]
  
  # IC of neighbouring trees at each direction
  # normalize to make all coefficients-squared add up to one
  BIC <- sapply(seq_along(dirs),
                function(dir) ifelse(is.na(Bneigh[, dir]), NA, wdir.dist[dir]))
  BIC <- BIC / sqrt(apply(BIC**2, 1, sum, na.rm = TRUE))
  
  # check
  stopifnot(all(sapply(apply(BIC**2, 1, sum, na.rm = TRUE), all.equal, 1)))
  
  # cbind the Intensity of Competition and the neighbours idx
  B16 <- cbind(BIC, Bneigh)
  B16[is.na(B16)] <- 0
  
  inc.mat <- matrix.short16(B16)
  
  ## Build the spatial object and add the 'competition' class
  spatial.call <- mc
  spatial.call[[1]] <- as.symbol('spatial')
  spatial.call <- spatial.call[-grep('decay|autofill', names(mc))]
  spatial.call$incidence <- inc.mat
  ans <- eval(spatial.call, parent.frame())
  class(ans) <- c('competition', class(ans))
  
  return(ans)
}


#' Build a permanent-environmental competition model
#' 
#' Given the coordinates of observations, and the competition decay parameter, 
#' build a \code{permanent_environmental_competition} model.
#' 
#' @inheritParams build_grid
#' @inheritParams competition
#'   
#' @return An object inheriting from \link{\code{competition}} with the
#'   incidence and structure matrices for the random effect.
#' @examples 
#' dat <- data.frame(id   = 1:5,
#'                   sire = c(11, 11, 2, 3, 2),
#'                   dam  = c(12, NA, 1, 12, 1),
#'                   x    = c(rep(1:2, times = 2), 3),
#'                   y    = c(rep(1:2, each = 2), 3))
#' permanent_environmental_competition(coord = dat[, c('x', 'y')], decay = 2)
permanent_environmental_competition <- function(coordinates,
                                                decay,
                                                autofill = TRUE) {
  
  ## Checks
  
  ## Use a diagonal covariance for pec
  cov.mat <- diag(nrow(coordinates))

  ## Build competition object  
  ans <- competition(coordinates = coordinates, 
                     covariance  = cov.mat,
                     decay       = decay,
                     autofill    = autofill)
  
  ## Set class
  class(ans) <- c('permanent_environmental_competition',
                  class(ans))
  
  return(ans)
}

