#' Build competition model
#' 
build.competition.model <- function (genetic) {
  
  # Checks
  
  # Original coordinates
  coord0 <- as.data.frame(sapply(coord, as.numeric))
    
  # lattice of spatial locations
  # possibly with automatic filling of empty rows or columns
  pos <- loc_grid(coord, autofill)
  
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
  Uinv <- kronecker(Q1d[[2]], Q1d[[1]])
  dimUinv <- dim(Uinv)
  stopifnot(identical(dimUinv[1], dimUinv[2]))
  
  # Map data coordinates with corresponding index of the Q matrix
  matrix2vec <- function(x, nx = pos.length[1], ny = pos.length[2]) {
    map <- matrix(1:(nx*ny), nx, ny)
    return(apply(x, 1, function(y) map[y[1], y[2]]))
  }
  data.ordering <- matrix2vec(sapply(coord, as.integer))
  
  
  # Incidence matrix
  # It is possible that the pedigree has been recoded/reordered
  # In that case, we need to recode the data file id codes as well
  if( !is.null(attr(genetic$pedigree, 'map')) )
    idx <- attr(genetic$pedigree, 'map')[genetic$id]
  else
    idx <- genetic$id
  
  B = idx
  if( genetic$model == 'competition') {
    # put the Intensity of Competition here
    B <- cbind(B, IoC)
    # put the neighbourhood structure here
    # Each individual has (up to) eight neighbours.
    # Moving clockwise from North, this matrix stores the id
    # of the corresponding neighbourg, or zeroes if missing
    Bnest 
  }
  
    
  # Scaling so that the characteristic marginal variance equals 1/sigma^2
  # Sorbye and Rue (2014)
  # Caveat: I need to invert the matrix here
  # Is there a way of finding the characteristic marginal variance
  # from the precision matrix?
  B <- sparseMatrix(i = 1:length(data.ordering),
                    j = data.ordering,
                    x = 1)
  
  Uinv <- Uinv * gmean(Matrix::diag(B %*% solve(Uinv) %*% Matrix::t(B)))
  
  # Store only the lower triangle of the symmetric matrix Q in triplet form
  trilUinv <- as(tril(Uinv), 'dgTMatrix')
  
  # Coordinates for the full grid
  #   coord.1d <- function(nx, levels) {
  #     as.numeric(factor(1:nx, labels = levels))
  #     # Isn't this the same as 1:nx?
  #   }
  #   plot.grid <- expand.grid(mapply(coord.1d, pos.length, lapply(coord, levels),
  #                                   SIMPLIFY = FALSE))
  plot.grid <- expand.grid(pos)
  plotting <- list(grid = plot.grid)
  return(list(param = genetic$competition_decay,
              coord = coord0,
              # Either the id
              B = data.ordering,  
              U = cbind(trilUinv@i + 1, trilUinv@j + 1, trilUinv@x),
              Utype = 'precision',
              plotting = plotting))
}



# 'move' an arrangement in a given direction
neighbours.at <- function(x, dir) UseMethod('neighbours.at')
neighbours.at.matrix <-function(x, dir) {
  stopifnot(is.character(dir))
  if( length(dir) == 1 ){
    dimx <- dim(x)
    switch(dir,
           N = rbind(NA, x[-dimx[1], ]),
           S = rbind(x[-1, ], NA),
           E = cbind(x[, -1], NA),
           W = cbind(NA, x[, -dimx[2]]),
           NE = cbind(rbind(NA, x[-dimx[1], -1]), NA),
           SE = cbind(rbind(x[-1, -1], NA), NA),
           SW = cbind(NA, rbind(x[-1, -dimx[2]], NA)),
           NW = cbind(NA, rbind(NA, x[-dimx[1], -dimx[2]]))
    )
  } else {
    sapply(dir, function(d) neighbours.at(x, d))
  }
}
neighbours.at.list <- function(x, dir) {
  # Check that it is a list of matrices
  stopifnot(all(lapply(x, class) == 'matrix'))
  lapply(x, neighbours.at, dir)
}

