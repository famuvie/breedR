#' Build additive direct and competition models
#' 
build.genetic.model <- function (genetic) {
  
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
    coord0 <- as.data.frame(sapply(genetic$coord, as.integer))
    
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
    BIC <- sapply(1:length(dirs),
                  function(dir) ifelse(is.na(Bneigh[, dir]), NA, wdir.dist[dir]))
    BIC <- BIC / sqrt(apply(BIC**2, 1, sum, na.rm = TRUE))
    
    # check
    stopifnot(all(sapply(apply(BIC**2, 1, sum, na.rm = TRUE), all.equal, 1)))
    
    
      
    # cbind the Intensity of Competition and the neighbours idx
    B <- cbind(B, BIC, Bneigh)

    ans <- list(param = genetic$competition_decay,
                coord = coord0,
                B = B)
  }
  
  B[is.na(B)] <- 0
  return(c(ans, list(B = B)))
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

