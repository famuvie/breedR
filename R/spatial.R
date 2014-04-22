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