#' Build an autoregressive model
#' 
#' Given the coordinates of the observations, ...
build.ar.model <- function (coord, rho, autofill) {
  
  # Original coordinates
  coord0 <- as.data.frame(sapply(coord, as.numeric))
  
  # Sorted coordinates of rows and columns where there are
  # at least one observation (or missing)
  pos <- lapply(coord, function(x) sort(unique(as.numeric(x))))
  
  if( autofill ) {
    # Check wether there are full rows or columns missing
    pos <- mapply(fill_holes, pos, c('rows', 'columns'),
                  SIMPLIFY = FALSE)
  }
  # Otherwise, the positions are taken in order as if the 
  # spacing was regular
  
  # The coordinates as factors allows to find easily the number of 
  # different x and y positions, and the ordering
  coord <- as.data.frame(mapply(factor, as.data.frame(coord), pos, SIMPLIFY = FALSE))
  
  # Number of different locations in rows and cols
  pos.length <- sapply(pos, length)
  
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
  data.ordering <- matrix2vec(sapply(coord, as.integer))
  
  
  # Coordinates for the full grid
  #   coord.1d <- function(nx, levels) {
  #     as.numeric(factor(1:nx, labels = levels))
  #     # Isn't this the same as 1:nx?
  #   }
  #   plot.grid <- expand.grid(mapply(coord.1d, pos.length, lapply(coord, levels),
  #                                   SIMPLIFY = FALSE))
  plot.grid <- expand.grid(pos)
  plotting <- list(grid = plot.grid)
  return(list(param = rho,
              coord = coord0,
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
  rho <- as.data.frame(rho)
  # If this function was called, at least one of the parameters is NA (unknown)
  stopifnot( length(rho) == 2)
  # We start with a very rough approximation
  rho.values <- c(-8, -2, 2, 8)/10
  
  set.values <- function(r) if(all(is.na(r))) rho.values else r
  grid <- expand.grid(lapply(rho, set.values))
  names(grid) <- c('rho_r', 'rho_c')
  return(grid)
}

#' Fill a hole in a given position
fill_hole <- function(x, idx, n, sep, label) {
  filling <- x[idx] + sep*(1:n)
  if( !all.equal(x[idx + 1], x[idx] + sep*(n+1)) )
    stop(paste('There is a hole in the', label,
               'which is not a multiple of the',
               'separation between individuals.\n',
               'Please introduce missing observations',
               'as needed and use autofill = FALSE\n'))
  return(c(head(x, idx), filling, tail(x, -idx)))
}


#' Find and fill all the holes in a vector
fill_holes <- function(x, label) {
  
  # Distance between lines of individuals
  dif <- diff(x)
  
  # Check whether there is a "regular spacing"
  # defined as the spacing between at least 80% of the individuals
  if( diff(quantile(dif, probs = c(.1, .9))) != 0 )
    stop("This does not seem to be a regular grid.\n",
         "The spacing between", label, "should be the same for",
         "at least the 80% of the cases.\n",
         "You can override this check with autofill = FALSE.\n")
  
  # Otherwise, the grid spacing can be defined as the median sepparation
  sep <- median(dif)
  
  # Check whether there are some "holes"
  if( min(dif) != max(dif) ) {
    holes <- which(dif > sep)
    # The hole can be either larger or shorter than the standard sep
    hole_sizes <- dif[holes] - 1
    
    for( i in 1:length(holes) ) {
      x <- fill_hole(x, holes[i], hole_sizes[i] , sep, label)
      holes[i+1] <- holes[i+1] + hole_sizes[i]
    }
  }
  
  return(x)
}

