#' @describeIn get_structure Return the structure matrices of all structured random effects
#' @export
get_structure.breedR <- function (x) {
  
  ## Only works for refactored effects
  eg.idx <- vapply(x$effects, inherits, TRUE, 'effect_group')
  
  sml <- lapply(x$effects[eg.idx], get_structure.effect_group)
  
  return(sml)
}


#' @describeIn get_structure Check that all elements share the same structure and return it.
#' @export 
get_structure.effect_group <- function(x) {
  
  ## get the incidence matrices for all the subeffects
  str.list <- lapply(x$effects, get_structure.breedr_effect)

  ## confirm they have all the same structure
  if (nelem <- length(str.list) > 1) {
    str.types <- vapply(str.list, attr, '','type')
    
    ## Structure types can be either covariance or precision matrices
    legal.types <- c('covariance', 'precision')
    stopifnot(all(unique(str.types) %in% legal.types))
        
    if (length(unique(str.types)) == 1) {
      ## Compare matrices all of the same type
      stopifnot(all.equal(str.list[1], str.list[-1]))
      stopifnot(isTRUE(all.equal(str.list[[1]], str.list[[-1]])))
    } else {
      ## Compare matrices of the corresponding types
      str.list.cov <- str.list[str.types == 'covariance']
      str.list.prec <- str.list[str.types == 'precision']

      if (length(str.list.cov) > 1)
        stopifnot(all.equalx(str.list.cov[1], str.list.cov[-1]))
      if (length(str.list.prec) > 1)
        stopifnot(all.equalx(str.list.prec[1], str.list.prec[-1]))
      
      ## Compare one covariance with one inverted precision
      ## Converting to standard matrix format, as solving often 
      ## results in a dense matrix
      str.cov1 <- as(str.list.cov[[1]], 'matrix')
      str.cov2 <- as(Matrix::solve(str.list.prec[[1]]), 'matrix')
      stopifnot(all.equal(str.cov1, str.cov2))
    }
  }
  
  ## Return the structure of the first element as a representative
  return(str.list[[1]])
}



#' @describeIn get_structure Return the structure matrix with an attribute
#'   indicating its \code{type}.
#' @export
get_structure.breedr_effect <- function(x) {
  ans <- structure(x$structure.matrix,
                   type = x$structure.type)
  return(ans)
}
