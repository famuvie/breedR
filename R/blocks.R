#' Build a blocks model
#' 
#' Given the coordinates of the observations, the *factor* identifying blocks,
#' and the logical autofill, build the incidence and covariance matrices of a
#' blocks model.
#' 
#' @param coord matrix(-like) of observation coordinates
#' @param id factor of the same length as observations, giving the block id for
#'   each observation.
#' @param autofill logical. If \code{TRUE} (default) it will try to fill gaps in
#'   the rows or columns. Otherwise, it will treat gaps the same way as adjacent
#'   rows or columns.
breedR_blocks <- function (coord,
                           id,
                           autofill = TRUE) {
  
  ## Checks
  if (!is.factor(id))  id <- as.factor(id)
  
  # Consider matrix-like coordinates
  coordinates <- as.data.frame(coord)
  
  ## Encompassing grid
  grid <- build_grid(coordinates, autofill)
  
  # How to "fill-in" the missing locations with the right block number?
  # Not trivial. Consider the most common level among the neighbors?
  # For the moment, don't fill anything.
  
  # Structure matrix for the blocks (identity)
  # (needed by vcov())
  n.blocks <- nlevels(id)
  cov.mat <- diag(n.blocks)
  
  inc.mat <- Matrix::sparseMatrix(i = seq.int(id),
                                  j = as.numeric(id),
                                  x = 1)

  ## Build the spatial effect, return the autoregressive parameters
  ## and further specify the ar class
  ans <- spatial(coordinates, incidence = inc.mat, covariance = cov.mat)
  ans$param <- list(n.blocks = n.blocks)
  attr(ans, 'grid') <- grid
  class(ans) <- c('blocks', class(ans))
  
  return(ans)
  
}

