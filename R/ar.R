#' Build an autoregressive model
#' 
#' Given the coordinates of the observations, the autocorrelation parameters
#' and the autofill logical value, computes the incidence 
#' matrix B and the covariance matrix U
#' 
#' @param coord matrix(-like) of observation coordinates
#' @param rho numeric. Vector of length two with the autocorrelation parameter 
#'   in each dimension, with values in the open interval (-1, 1).
#' @param autofill logical. If \code{TRUE} (default) it will try to fill gaps in
#'   the rows or columns. Otherwise, it will treat gaps the same way as adjacent
#'   rows or columns.
#'   
#' @return a list with
#' - the autocorrelation parameters
#' - the coordinates original coordinates as a data frame
#' - the incidence matrix (encoded as a vector)
#' - the covariance matrix (of the full grid, in triplet form)
breedr_ar <- function (coord,
                       rho,
                       autofill = TRUE) {
  
  # Checks
  if( !all(abs(rho) < 1) )
    stop('The autoregressive parameters rho must be strictly in (-1, 1).\n')

  # Consider matrix-like coordinates
  coordinates <- as.data.frame(coord)
  
  ## Encompassing grid
  grid <- build_grid(coordinates, autofill)
  
  # Check for regular grid
  if (!grid$regular)
    stop('The AR model can only be fitted to regular grids.')
  
  # Check: the grid should be actually 2d
  if( !all(grid$length > 2) )
    stop('Are you kidding? This is a line!')
  
  ## Incidence matrix
  inc.mat <- Matrix::sparseMatrix(i = seq_along(grid$idx),
                                  j = grid$idx,
                                  x = 1,
                                  dims = c(length(grid$idx),
                                           prod(grid$length)))
  
  
  # Precision matrix for the AR1(rho_x) x AR1(rho_y) process
  # when locations are stacked following the standards of R:
  # by columns: first vary x and then y
  # (1, 1), (2, 1), ..., (n_x, 1), (1, 2), ..., (n_x, 2), ...
  # Only the lower triangle
  Q1d <- mapply(build.AR1d, grid$length, rho)
  Uinv <- kronecker(Q1d[[2]], Q1d[[1]])
  dimUinv <- dim(Uinv)
  stopifnot(identical(dimUinv[1], dimUinv[2]))
  
  
  # Scaling so that the characteristic marginal variance equals 1/sigma^2
  # Sorbye and Rue (2014)
  # In principle I would need to invert the matrix here
  # But we can derive the characteristic marginal variance
  # from the precision matrix as follows:
  # In the AR1xAR1 model, B is a permutation matrix, thus 
  # diag(B·U·B') = diag(U)
  # On the other hand, U = Q2inv %x% Q1inv, the kronecker product of the
  # one-dimentional autoregressive *covariance* matrices, which have constant
  # diagonal equal to 1/(1-rho^2).
  # An element in the diagonal of the kronecker product is the scalar product
  # of the corresponding elements in the diagonals (which are constant)
  # So, diag(U) is constant, and the geometric mean turns out to be the 
  # product of 1/(1-rho**2) for both rhos.
  scaling <- prod(1/(1-rho**2))
  Uinv <- Uinv * scaling
  
  ## Build the spatial effect, return the autoregressive parameters
  ## and further specify the ar class
  ans <- spatial(coordinates, incidence = inc.mat, precision = Uinv)
  ans$param <- list(rho = rho)
  attr(ans, 'grid') <- grid
  class(ans) <- c('ar', class(ans))
  
  return(ans)
}


# Build an evaluation grid for the autoregressive parameters of rows and columns
# 
# One or both autoregressive parameters may be unknown.
# In that case we need to fit the model with several values of rho_r and rho_c,
# in order to estimate the most likely values.
# This functions provides an evaluation grid of values.
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

# Build precision matrix of AR1(rho) in the line
build.AR1d <- function(n, x) {
  temp <- diag(c(1, rep(1 + x^2, n-2), 1))
  subdiag <- rbind(0, cbind(diag(-x, n-1), 0))
  return(as(Matrix::Matrix(temp + subdiag + t(subdiag), sparse = TRUE), 'dgTMatrix'))
}
