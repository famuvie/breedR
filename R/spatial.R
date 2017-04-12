#' Build a spatial model
#' 
#' Check conformity of arguments and return a \code{spatial} object.
#' 
#' @inheritParams random
#' @inheritParams build_grid
#' @return A list with elements \code{coordinates}, \code{incidence.matrix},
#'   \code{structure.matrix} and \code{structure.type}, which is a string
#'   indicating either \code{covariance} or \code{precision}.
spatial <- function(coordinates, incidence, covariance, precision) {
  
  mc <- match.call()
  arg.list <- as.list(mc)[-1]
  
  ## check conformance
  if (nrow(coordinates) != nrow(incidence))
    stop('The incidence matrix should have as many rows as observations.')
  
  ## Build the random effect, and further specify the spatial class
  random.args <- lapply(arg.list[-1], eval, parent.frame())
  ans <- do.call('random', random.args)
  ans$coordinates <- coordinates
  class(ans) <- c('spatial', class(ans))

  return(ans)
}



#' Lattice of spatial locations
#' 
#' Returns a list row and column coordinates of observations
#' 
#' @inheritParams build_grid
#' @return list of row and column coordinates of spatial units
#'   
#'   This functions converts observational coordinates into spatial coordinates.
#'   While in a observational dataset there can possibly be many or missing 
#'   observations at one single location, this function returns the list of row
#'   and column coordinates of a grid which contains all the observations
loc_grid <- function (coordinates, autofill) {
  # Sorted coordinates of rows and columns where there are
  # at least one observation (or missing)
  pos <- lapply(as.data.frame(coordinates),
                function(x) sort(unique(as.numeric(x))))
  
  if( autofill ) {
    # Check wether there are full rows or columns missing
    pos <- mapply(fill_holes, pos, c('rows', 'columns'),
                  SIMPLIFY = FALSE)
  }
  # Otherwise, the positions are taken in order as if the 
  # spacing was regular
  
  return(pos)
}



#' Find and fill all the holes in a vector
#' 
#' @param x numeric. Vector of increasing coordinates.
#' @param label character. A name like 'rows' or 'x'.
#' @importFrom stats quantile median
fill_holes <- function(x, label) {
  
  if(length(x) > 1){
    # Distance between lines of individuals
    dif <- diff(x)
    
    # Check whether there is a "regular spacing"
    # defined as the spacing between at least 60% of the individuals
    if( !isTRUE(all.equal(diff(quantile(dif, probs = c(.1, .6))), 0,
                          check.attributes = FALSE)) )
      stop("This does not seem to be a regular grid.\n",
           "The spacing between ", label, " should be the same for ",
           "at least the 60% of the cases.\n",
           "You can override this check with autofill = FALSE.\n")
    
    # Otherwise, the grid spacing can be defined as the median sepparation
    sep <- median(dif)
    
    # Check whether there are some "holes"
    if( !isTRUE(all.equal(min(dif), max(dif))) ) {
      holes <- which(dif > sep)
      # The hole can be either larger or shorter than the standard sep
      hole_sizes <- dif[holes]/sep - 1
      
      for( i in seq_along(holes) ) {
        x <- fill_hole(x, holes[i], hole_sizes[i] , sep, label)
        holes <- holes + hole_sizes[i]
      }
    }
  }
  
  return(x)
}


# Fill a hole in a given position
fill_hole <- function(x, idx, n, sep, label) {
  filling <- x[idx] + sep*(1:n)
  if( !isTRUE(all.equal(x[idx + 1], x[idx] + sep*(n+1))) )
    stop(paste('There is a hole in the', label,
               'which is not a multiple of the',
               'separation between individuals.\n',
               'Please introduce missing observations',
               'as needed and use autofill = FALSE\n'))
  return(c(head(x, idx), filling, tail(x, -idx)))
}


#' Build an encompassing grid
#' 
#' Build the minimal regularly-spaced grid containing a given set of points.
#' 
#' Note that \code{autofill = FALSE} virtually removes the empty lines,
#' considering the spacing as constant.
#' 
#' @param coordinates two-column matrix-like set of row and column coordinates 
#'   of observational units
#' @param autofill logical. If TRUE (default) it will try to fill missing rows
#'   or columns with missing observations. Otherwise, will treat individuals as
#'   neighbours even if they are across an empty line.
#'   
#' @return The parameters defining the grid, and the index of the observed 
#'   coordinates in the grid seen as a vector. More specifically, \describe{ 
#'   \item{origin}{the coordinates of the \emph{first} (with smallest row and 
#'   column values) point in the grid} \item{step}{the separation between rows 
#'   and columns} \item{length}{the number of points in each dimension} 
#'   \item{idx}{the index of each observation in the vectorized grid} }
build_grid <- function (coordinates, autofill = TRUE) {
  
  stopifnot(ncol(coordinates) == 2)
  
  ## Original coordinates
  coord0 <- as.data.frame(coordinates)
  
  ## spatial locations by rows and columns
  ## possibly with automatic filling of empty rows or columns
  pos <- loc_grid(coord0, autofill)
  
  # The coordinates as factors allows to find easily the number of 
  # different x and y positions, and the ordering
  coord <- as.data.frame(mapply(factor,
                                coord0,
                                pos,
                                SIMPLIFY = FALSE))
  
  # Number of different locations in rows and cols
  pos.length <- vapply(pos, length, 1)
  
  # Spacing between trees in rows and columns
  # If autofill, then diff(x) is constant and equal to its min
  # If !autofill, then there might be some holes and the spacing
  # should be defined as the minimal separation.
  pos.step <- vapply(pos, function(x) min(diff(x)), 1)
  
  # Map data coordinates with corresponding index of the Q matrix
  matrix2vec <- function(x, nx = pos.length[1], ny = pos.length[2]) {
    map <- matrix(1:(nx*ny), nx, ny)
    return(apply(x, 1, function(y) map[y[1], y[2]]))
  }
  ord <- matrix2vec(sapply(coord, as.integer))
  
  ## Check for regular grid
  # if regular, n_x \times n_y ~ n_obs
  # if irregular, n_x \times n_y ~ n_obs^2
  # We assume it is regular if
  # n_x \times n_y < n_obs + (n_obs^2 - n_obs)/4
  looks.irregular <- prod(pos.length) > nrow(coordinates)*(nrow(coordinates) + 3)/4
  
  return(list(origin   = vapply(pos, head, 1, 1),
              step     = pos.step,
              length   = pos.length,
              idx      = ord,
              regular  = !looks.irregular,
              autofill = autofill))
}