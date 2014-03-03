#' Build an autoregressive model
#' 
#' Given the coordinates of the observations, ...
build.ar.model <- function (coord, rho) {
  # The coordinates as factors allows to find easily the number of 
  # different x and y positions, and the ordering
  pos <- lapply(as.data.frame(coord), factor)
  
  # Number of different locations in rows and cols
  pos.length <- sapply(pos, function(x) length(levels(x)))
  
  # Check: is this a regular grid?
  # if regular, n_x \times n_y ~ n_obs
  # if irregular, n_x \times n_y ~ n_obs^2
  # We assume it is regular if
  # n_x \times n_y < n_obs + (n_obs^2 - n_obs)/4
  if(prod(pos.length) > nrow(coord)*(nrow(coord) + 3)/4) 
    stop('The AR model can only be fitted to regular grids.')
  
  # Check: the grid should be actually 2d
  if( !all(pos.length > 2) )
    stop('Are you kidding? This is a line!')
  
  # Build precision matrix of AR1(rho) in the line
  build.AR1d <- function(n, x) {
    temp <- diag(c(1, rep(1 + x^2, n-2), 1))
    subdiag <- rbind(0, cbind(diag(-x, n-1), 0))
    return(as(Matrix(temp + subdiag + t(subdiag), sparse = TRUE), 'dgTMatrix'))
  }
  Q1d <- mapply(build.AR1d, pos.length, rho)

  # Precision matrix for the AR1(rho_x) x AR1(rho_y) process
  # when locations are stacked following the standards of R:
  # by columns: first vary x and then y
  # (1, 1), (2, 1), ..., (n_x, 1), (1, 2), ..., (n_x, 2), ...
  # Only the lower triangle
  Q <- tril(kronecker(Q1d[[2]], Q1d[[1]]))
  # Map data coordinates with corresponding index of the Q matrix
  matrix2vec <- function(x, nx = pos.length[1], ny = pos.length[2]) {
    map <- matrix(1:(nx*ny), nx, ny)
    return(apply(x, 1, function(y) map[y[1], y[2]]))
  }
  data.ordering <- matrix2vec(sapply(pos, as.integer))
  
  
  # Coordinates for the full grid
  coord.1d <- function(nx, levels) {
    as.numeric(factor(1:nx, labels = levels))
  }
  plot.grid <- expand.grid(mapply(coord.1d, pos.length, lapply(pos, levels)))
  
  plotting <- list(grid = plot.grid)
  return(list(coord = coord,
              B = data.ordering,
              U = cbind(Q@i + 1, Q@j + 1, Q@x),
              plotting = plotting))
}

#' Build an evaluation grid for the autoregressive parameters of rows and columns
#' 
#' One or both autoregressive parameters may be unknown.
#' In that case we need to fit the model with several values of rho_r and rho_c,
#' in order to estimate the most likely values.
#' This functions provides an evaluation grid of values.
build.AR.rho.grid <- function(rho) {
  # If this function was called, at least one of the parameters is NA (unknown)
  stopifnot( length(rho) == 2)
  # We start with a very rough approximation
  rho.values <- seq(-.5, .5, length = 2)
  
  set.values <- function(r) if(is.na(r)) rho.values else r
  grid <- expand.grid(lapply(rho, set.values))
  names(grid) <- c('rho_r', 'rho_c')
  return(grid)
}