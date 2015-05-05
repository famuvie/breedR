#' Build a spatial model
#' 
#' Check conformity of arguments and return a \code{spatial} object.
#' 
#' @param coord two-column matrix-like object with observation coordinates
#' @inheritParams random
#' @return A list with elements \code{coordinates}, \code{incidence.matrix},
#'   \code{structure.matrix} and \code{structure.type}, which is a string
#'   indicating either \code{covariance} or \code{precision}.
spatial <- function(coord, incidence, covariance, precision) {
  
  mc <- match.call()
  arg.list <- as.list(mc)[-1]
  
  ## check conformance
  if (nrow(coord) != nrow(incidence))
    stop('The incidence matrix should have as many rows as observations.')
  
  ## Build the random effect, and further specify the spatial class
  random.args <- lapply(arg.list[-1], eval, parent.frame())
  ans <- do.call('random', random.args)
  ans$coordinates <- coord
  class(ans) <- c('spatial', class(ans))

  return(ans)
}



#' Lattice of spatial locations
#' 
#' Returns a list row and column coordinates of observations
#' 
#' @param coord matrix, data.frame or list of row and column coordinates of 
#'   observational units
#' @param autofill locigcal, whether to perform automatic filling of empty 
#'   rows/columns
#' @return list of row and column coordinates of spatial units
#'   
#'   This functions converts observational coordinates into spatial coordinates.
#'   While in a observational dataset there can possibly be many or missing 
#'   observations at one single location, this function returns the list of row
#'   and column coordinates of a grid which contains all the observations
loc_grid <- function (coord, autofill) {
  # Sorted coordinates of rows and columns where there are
  # at least one observation (or missing)
  pos <- lapply(as.data.frame(coord),
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
fill_holes <- function(x, label) {
  
  if(length(x) > 1){
    # Distance between lines of individuals
    dif <- diff(x)
    
    # Check whether there is a "regular spacing"
    # defined as the spacing between at least 60% of the individuals
    if( !isTRUE(all.equal(diff(quantile(dif, probs = c(.1, .6))), 0, check.attributes = FALSE)) )
      stop("This does not seem to be a regular grid.\n",
           "The spacing between", label, "should be the same for ",
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


#' Fill a hole in a given position
#' 
#' @param x numeric vector
#' @param idx index of last observation just before a hole
#' @param n size of hole in terms of number of missing observations
#' @param sep spacing between observations
#' @param label character. An identificative name.
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
