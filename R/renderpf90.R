#' @describeIn renderpf90 For unknown classes, just returns the object untouched
#' @export
renderpf90.default <- function(x) x


#' Render a sparse matrix into non-zero values and column indices
#' 
#' This function provides a representation of a sparse matrix suited to 
#' progsf90.
#' 
#' For each row, it keeps the non-zero elements and their respective column 
#' index. The gaps are filled with zeros.
#' 
#' @return A matrix with a number of colums equal to twice the maximum number of
#'   non-zero elements in one row. The first half of the columns are the 
#'   non-zero values (except for filling-in) while the second half are the 
#'   column indices.
#' @inheritParams renderpf90
#' @family renderpf90
#' @export
renderpf90.matrix <- function(x) {
  
  which_n_where <- function(x, n) {
    w <- which(x > 0)
    x <- x[w]
    if (sup <- n - length(w)) {
      w <- c(w, rep(0, sup))
      x <- c(x, rep(0, sup))
    }
    
    c(x, w)
  }
  
  x.nz <- x > 0
  ## Effective number of columns
  ncol_eff <- max(apply(x.nz, 1, sum))
  t(apply(x, 1, which_n_where, ncol_eff))
}


#' @describeIn renderpf90 Render a full \code{breedr_modelframe}
#' @export
renderpf90.breedr_modelframe <- function(x) {
  
  ## Until the refactoring is completed, not all effects in the list
  ## will be breedr_effects or effect_groups
  ## Those which are not remain untouched
  
  xpf90 <- lapply(x, renderpf90)
  
  ## The dimension of the response will be translated here from progsf90
  ## This determines the initial position of the effects in the data file

  ## Positions of the 'virtual' effects 
  ## (to appear in the EFFECTS section in progsf90)
  ## Temporarily, only needed for refactored effect_groups
  aux.idx <- which(sapply(x, inherits, 'effect_group'))
  offset <- 0
  for (i in aux.idx) {
    pos.length <- length(xpf90[[i]]$levels)
    stopifnot(all.equal(pos.length, length(xpf90[[i]]$nest)))
    
    xpf90[[i]]$pos <- x[[i]]$pos.head - 1 + seq.int(from = 1, to = pos.length) + offset
    xpf90[[i]]$type <- paste(xpf90[[i]]$type, x[[i]]$pos.head - 1 + xpf90[[i]]$nest + offset)
    ## Make sure file_name is unique
    xpf90[[i]]$file_name <- paste(xpf90[[i]]$file_name, names(xpf90)[i], sep = '_')
    
    offset <- offset + 2*pos.length
    
    
  }
  
  return(xpf90)
}
  

#' @describeIn renderpf90 Render groups of effects into pf90 code
#' @export
renderpf90.effect_group <- function(x) {

  ## Render each effect individually
  
  #   if (length(x$effects) == 1) {  # Unnecessary (although faster)
  #     ans <- renderpf90(x$effects[[1]])
  #     ans$var <- x$cov.ini
  #     return(ans)
  #   }
  aux <- lapply(x$effects, renderpf90)
  
  ## Check that the model name and file are unique,
  ## All elements in the group must have the same structure
  concatenate <- function(x, item) 
    do.call('c', lapply(aux, function(x) x[[item]]))
  
  models <- unique(concatenate(aux, 'model'))
  fnames <- unique(concatenate(aux, 'file_name'))
  if (length(models) != 1 | length(fnames) != 1)
    stop('progsf90 only admits random groups with the same structure.')
  
  ## TODO
  ## Assume that all structure matrices (in file) are equal
  ## Need to check for that.
  
  ## DATA file
  ## incidence matrix of the group rendered for pf90
  dat <- renderpf90.matrix(model.matrix.effect_group(x))
  
  
  ## Concatenate parameters and bind data column-wise
  ## as this is a restriction in progsf90
  list(
    levels    = concatenate(aux, 'levels'),
    type      = concatenate(aux, 'type'),
    nest      = concatenate(aux, 'nest'),
    model     = models,
    file_name = fnames,
    file      = aux[[1]]$file,  # Assuming all are identical
    var       = x$cov.ini,
    data      = dat)
}



#' @details For the \code{generic} method, all matrices are converted to plain 
#'   matrix-class, for exporting to files.
#' @return For the \code{generic} method, number of levels and type for each
#'   'virtual' effect; model \code{user_file} or \code{user_file_i} as
#'   appropriate; a file name and its content.
#' @describeIn renderpf90 Compute the parameters of a progsf90 representation of
#'   a generic effect.
#' @export
renderpf90.generic <- function(x) {
  
  ## Incidence matrix in 'reduced' form
  mm <- as.matrix(model.matrix(x))
  mmpf90 <- renderpf90.matrix(mm)
  
  ## Number of 'virtual' effects
  ## i.e. half the number of columns in the reduced matrix
  nve <- ncol(mmpf90)/2
  
  ## Structure matrix (either covariance or precision)
  sm <- vcov(x)
  smpf90 <- as.triplet(sm)
  
  ## Number of levels for each of the effects in progsf90
  ans <- list(
    levels    = c(rep(0, nve - 1), nrow(sm)),
    type      = rep('cov', nve),
    nest      = nve + 1:nve,
    model     = ifelse(attr(sm, 'inverse'), 'user_file', 'user_file_i'),
    file_name = 'generic',
    file      = smpf90
  )
  return(ans)
}

#' Represent a symmetric matrix in triplet format
as.triplet <- function(x) {
  xsp <- Matrix::tril(as(as.matrix(x), 'dgTMatrix'))
  # Note: The Matrix package counts rows and columns starting from zero
  # Thus, I add 1 to the corresponding columns
  cbind(xsp@i + 1, xsp@j + 1, xsp@x)
}