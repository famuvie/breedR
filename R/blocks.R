# Build a blocks model
# 
# Given the coordinates of the observations,
# the *factor* identifying blocks, and the logical autofill
build.blocks.model <- function (coord, id, autofill) {
  
  stopifnot(is.factor(id))  # This should not happen
  
  # Original coordinates
  coord0 <- as.data.frame(sapply(coord, as.numeric))
    
  # lattice of spatial locations
  # possibly with automatic filling of empty rows or columns
  pos <- loc_grid(coord, autofill)
  
  # The coordinates as factors allows to find easily the number of 
  # different x and y positions, and the ordering
  coord <- as.data.frame(mapply(factor,
                                as.data.frame(coord),
                                pos, SIMPLIFY = FALSE))
  
  # Number of different locations in rows and cols
  pos.length <- sapply(pos, length)
  
  # This will map the observed positions in the full grid
  matrix2vec <- function(x, nx = pos.length[1], ny = pos.length[2]) {
    map <- matrix(1:(nx*ny), nx, ny)
    return(apply(x, 1, function(y) map[y[1], y[2]]))
  }
  data.ordering <- matrix2vec(sapply(coord, as.integer))
  # How to "fill-in" the missing locations with the right block number?
  # Not trivial. Look the most common level among the neighbors?
  # For the moment, don't fill anything.
  
  # Structure matrix for the blocks (identity)
  # (needed by vcov())
  n.blocks <- nlevels(id)
  U <- cbind(1:n.blocks, 1:n.blocks, 1)
  
  plot.grid <- expand.grid(pos)
  inc.mat <- Matrix::sparseMatrix(i = data.ordering,
                                  j = as.numeric(id),
                                  x = 1)
  plotting <- list(grid = plot.grid,
                   inc.mat = inc.mat)
  return(list(coord = coord0,
              map = data.ordering,
              B = as.numeric(id),
              U = U,
              Utype = 'covariance',
              plotting = plotting))
}
