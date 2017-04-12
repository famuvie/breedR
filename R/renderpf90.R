#' @describeIn renderpf90 For unknown classes, just returns the object untouched
#'   with a warning.
#' @export
renderpf90.default <- function(x) {
  warning(paste('No method for',
                paste(class(x), collapse = ', '),
                'for renderpf90().'))
  return(x)
}


#' @describeIn renderpf90 For fixed effects.
#' @export
renderpf90.fixed <- function(x) {
  
  ## Do not use model.matrix() as this can be a factor
  mmx <- x$incidence.matrix
  
  if (is.factor(mmx)) {
    nl <- nlevels(mmx)
    type <- 'cross'
  } else {
    nl <- 1
    type <- 'cov'
  }
  
  return(list(pos = 1,
              levels = nl,
              type =type,
              nest = NA,
              data = as.matrix(as.numeric(mmx))))
}


#' @describeIn renderpf90 For diagonal effects. Assumed grouping variables.
#' @export
renderpf90.diagonal <- function(x) {
  
  ## Do not use model.matrix() as this can be a factor
  mmx <- x$incidence.matrix
  
  if (inherits(mmx, 'indMatrix')) {
    ## there are no missing values
    ## but there might be unobserved levels
    dat <- mmx@perm
  } else {
    ## general sparse matrix format: dgCMatrix
    stopifnot(all(mmx@x == 1))
    dp <- diff(mmx@p)
    j <- rep(seq_along(dp),dp)
    dat <- vector('integer', length = nrow(mmx))
    dat[mmx@i + 1] <- j
  }
  
  ans <- list(pos = 1,
              levels = ncol(mmx),
              type   = 'cross',
              nest = NA,
              model = 'diagonal',
              file = '',
              file_name = '',
              data = as.matrix(dat))

  return(ans)
}


#' Render a sparse matrix into non-zero values and column indices
#' 
#' This function provides a representation of a sparse matrix suited to 
#' progsf90.
#' 
#' For each row, it keeps the non-zero elements and their respective column 
#' index. The gaps are filled with zeros.
#' 
#' @return A matrix with a number of columns equal to twice the maximum number of
#'   non-zero elements in one row. The first half of the columns are the 
#'   non-zero values (except for filling-in) while the second half are the 
#'   column indices.
#'   
#'   For indicator matrices (i.e. each row has at most one non-zero value of 1),
#'   it returns a one-column matrix of the corresponding column indices.
#' @inheritParams renderpf90
#' @family renderpf90
#' @import Matrix
#' @export
renderpf90.matrix <- function(x) {
  
  ## Effective number of columns
  ## issue #37: excesive memory consumption
  ## using apply(x, 1, fun) where x is sparseMatrix implies
  ## a call to as.matrix(), potentially giving memory problems
  ncol_eff <- max(rowSums(x > 0))
  idx <- Matrix::which(x > 0, useNames = FALSE)  # non-zero locations

  ## non-zero matrix of positions and values sorted by position
  nzmat <- cbind(arrayInd(idx, .dim = dim(x)), x[idx])
  nzmat <- nzmat[ order(nzmat[, 1], nzmat[, 2]),]
  
  ## Fill-in zeros for spare columns
  fill_zero <- function(x, n) {
    if ((m <- length(x)) < n)
      x <- c(x, rep(0, n-m))
    return(x)
  }
  
  ## there might be some empty rows that
  ## we need to fill-in with zeros
  ans <- matrix(0, nrow = nrow(x), ncol = ncol_eff * 2)
  
  if (ncol_eff == 1) {
    ## Only one non-zero value per row
    stopifnot(all(!duplicated(nzmat[,1])))
    ans[nzmat[, 1], ] <- nzmat[, 3:2]
    
  } else {
    ## cbind non-zero values and their positions
    ## filling-in with zeros up to ncol_eff
    ans[unique(nzmat[, 1]), ] <- 
      cbind(do.call('rbind',
                    tapply(nzmat[, 3], nzmat[, 1], fill_zero, ncol_eff)),
            do.call('rbind',
                    tapply(nzmat[, 2], nzmat[, 1], fill_zero, ncol_eff)))
  }
  
  ## Exception: case of indicator matrix
  ## ncol_eff = 1, all values are 1 (or 0 if missing)
  ## return only column with positions
  if (ncol_eff == 1 && all(ans[, 1] <= 1))
    ans <- ans[, 2, drop = FALSE]
  
  return(ans)
}


#' @describeIn renderpf90 Render a full \code{breedr_modelframe}
#' @param ntraits integer. Number of traits in the model.
#' @param weights logical. Whether there is an additional column of weights.
#' @export
renderpf90.breedr_modelframe <- function(x, ntraits, weights) {
  
  xpf90 <- lapply(x, renderpf90)
  
  ## The dimension of the response (i.e. ntraits) will be translated here from
  ## progsf90 This determines the initial position of the effects in the data
  ## file, and the dimension of 'pos' and 'nest'.

  ## Make sure file_name are unique
  aux.idx <- which(sapply(x, inherits, 'effect_group'))
  for (i in aux.idx) {
    if ((fn <- xpf90[[i]]$file_name) != '') {
      xpf90[[i]]$file_name <- paste(fn,
                                    names(xpf90)[i], sep = '_')
    }
  }
  
  ## Positions of the 'virtual' effects in the data file
  ## (to appear in the EFFECTS section in progsf90)
  dat.widths <- sapply(xpf90, function(x) ncol(x$data))
  end.columns <- cumsum(c(response = ntraits + weights, dat.widths)) 
  offsets <- structure(head(end.columns, -1),
                       names = names(dat.widths))
  
  collapse_traits <- function(x, ntraits) {
    if (!is.na(x))
      paste(rep(x, ntraits), collapse = " ")
    else ''
  }

  for (i in seq_along(xpf90)) {
    xpf90[[i]]$pos <- vapply(offsets[[i]] + xpf90[[i]]$pos,
                          collapse_traits, "str", ntraits)
    xpf90[[i]]$nest <- vapply(offsets[[i]] + xpf90[[i]]$nest,
                           collapse_traits, "str", ntraits)
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
    unname(do.call('c', lapply(aux, function(x) x[[item]])))
  
  models <- unique(concatenate(aux, 'model'))
  fnames <- unique(concatenate(aux, 'file_name'))
  if (length(models) != 1 | length(fnames) != 1)
    stop('progsf90 only admits random groups with the same structure.')
  
  ## TODO
  ## Assume that all structure matrices (in file) are equal
  ## Need to check for that.
  
  ## DATA file
  ## cbind of rendered incidence matrices for each element of the group
  ## Note: it is not the same as the render of the group incidence matrix!
  dat.l <- lapply(lapply(x$effects, model.matrix), renderpf90.matrix)
  
  ## Update positions indices
  if ( length(aux) > 1) {
    dat.widths <- sapply(aux, function(x) ncol(x$data))
    end.columns <- cumsum(c(0, dat.widths)) 
    offsets <- structure(head(end.columns, -1),
                         names = names(dat.widths))
    for (i in seq_along(aux)) {
      aux[[i]]$pos <- offsets[[i]] + aux[[i]]$pos
      aux[[i]]$nest <- offsets[[i]] + aux[[i]]$nest
    }
  }
  
  ## Concatenate parameters and bind data column-wise
  ## as this is a restriction in progsf90
  list(
    pos       = concatenate(aux, 'pos'),
    levels    = concatenate(aux, 'levels'),
    type      = concatenate(aux, 'type'),
    nest      = concatenate(aux, 'nest'),
    model     = models,
    file_name = fnames,
    file      = aux[[1]]$file,  # Assuming all are identical
    var       = x$cov.ini,
    data      = do.call('cbind', dat.l))
}



#' @details For the \code{generic} class, all matrices are converted to plain 
#'   matrix-class, for exporting to files. The progsf90 model is either
#'   \code{user_file} or \code{user_file_i} depending on the type of structure
#'   matrix; i.e. respectively precision or covariance.
#' @describeIn renderpf90 Compute the parameters of a progsf90 representation of
#'   a generic effect.
#' @export
renderpf90.generic <- function(x) {
  
  ## Incidence matrix in 'reduced' form
  mmpf90 <- renderpf90.matrix(model.matrix(x))
  
  ## Number of 'virtual' effects
  ## i.e. half the number of columns in the reduced matrix
  ## except in the case of one column
  if ((nve <- nc <- ncol(mmpf90)) > 1)
    nve <- nc/2
  
  ## Structure matrix (either covariance or precision)
  sm <- vcov(x)
  smpf90 <- as.triplet(sm)
  
  ## Represent the data with a single factor
  ## if either the mm is one column (e.g. AR)
  ## or is a degenerate nested case: each row exactly one non-zero value
  ## (i.e. first column of 1). Note that there are two-column settings
  ## where some rows are full 0.
  single.column <- nc == 1 || (nc == 2 && all(mmpf90[, 1] == 1))
  
  if(single.column) {
    ## Reduced one-column form
    nest <- NA
    type <- 'cross'
    dat  <- mmpf90[, ncol(mmpf90), drop = FALSE] # last col (i.e. 1 or 2)
  } else {
    nest <- nve + 1:nve 
    type <- rep('cov', nve)
    dat  <- mmpf90
  }

  ## Number of levels for each of the effects in progsf90
  ans <- list(
    pos       = seq_len(nve),
    levels    = c(rep(0, nve - 1), nrow(sm)),
    type      = type,
    nest      = nest,
    model     = ifelse(attr(sm, 'inverse'), 'user_file', 'user_file_i'),
    file_name = 'generic',
    file      = smpf90,
    data      = dat
  )
  return(ans)
}

#' @describeIn renderpf90 Compute the parameters of a progsf90 representation of
#'   a additive_genetic_animal effect.
#' @export
renderpf90.additive_genetic_animal <- function(x) {
  
  ## pedigree in data.frame format
  ped.dat <- as.data.frame(get_pedigree(x))
  
  ans <- list(
    pos       = 1,
    levels    = nrow(ped.dat),
    type      = 'cross',
    nest      = NA,
    model     = 'add_animal', # both add_animal or competition
    file_name = 'pedigree',
    file      = ped.dat,
    data      = renderpf90.matrix(model.matrix(x))
    )

  return(ans)
}


#' @describeIn renderpf90 Compute the parameters of a progsf90 representation of
#'   a additive_genetic_competition effect.
#' @export
renderpf90.additive_genetic_competition <- function(x) {
  
  ## pedigree in data.frame format
  ped.dat <- as.data.frame(get_pedigree(x))
  
  ## maximum number of competitors
  n.comp <- 8
  
  ## dimension of the effect
  ef.dim <- nrow(ped.dat)
  
  ans <- list(
    pos       = seq_len(n.comp),
    levels    = c(rep(0, n.comp - 1), ef.dim),
    type      = rep('cov', n.comp),
    nest      = n.comp + seq_len(n.comp),
    model     = 'add_animal', # both add_animal or competition
    file_name = 'pedigree',
    file      = ped.dat,
    data      = renderpf90.matrix(model.matrix(x)))
  
  return(ans)
}


#' @describeIn renderpf90 Compute the parameters of a progsf90 representation of
#'   a permanent_environmental_competition effect. Has the same incidence matrix
#'   than the additive_genetic_competition effect but an unstructured covariance
#'   matrix.
#' @export
renderpf90.permanent_environmental_competition <- function(x) {
  
  ## incidence matrix
  mmx <- model.matrix(x)  
  
  ## maximum number of competitors
  n.comp <- 8
  
  ## dimension of the effect
  ef.dim <- ncol(mmx)
  
  ans <- list(
    pos       = seq_len(n.comp),
    levels    = c(rep(0, n.comp - 1), ef.dim),
    type      = rep('cov', n.comp),
    nest      = n.comp + seq_len(n.comp),
    model     = 'diagonal',
    file_name = '',
    file      = '',
    data      = renderpf90.matrix(model.matrix(x)))
  
  return(ans)
}


#' @details For the \code{splines} class, everything reduces to a generic effect
#'   with a covariance matrix
#' @describeIn renderpf90 Compute the parameters of a progsf90 representation of
#'   a splines effect.
#' @export
renderpf90.splines <- function(x) {
  
  ans <- renderpf90.generic(x)
  ans$file_name = 'splines'
  
  return(ans)
}



#' @details For the \code{blocks} class, everything reduces to a generic effect
#'   with a covariance matrix
#' @describeIn renderpf90 Compute the parameters of a progsf90 representation of
#'   a blocks effect.
#' @export
renderpf90.blocks <- function(x) {
  
  ## Do not use model.matrix() as this can be a factor
  mmx <- x$incidence.matrix
  
  ## general sparse matrix format: dgCMatrix
  stopifnot(all(mmx@x == 1))
  dp <- diff(mmx@p)
  j <- rep(seq_along(dp),dp)
  dat <- vector('integer', length = nrow(mmx))
  dat[mmx@i + 1] <- j
  
  ans <- list(pos = 1,
              levels = ncol(mmx),
              type   = 'cross',
              nest = NA,
              model = 'diagonal',
              file = '',
              file_name = '',
              data = as.matrix(dat))
  
  return(ans)
}


#' @details For the \code{ar} class, everything reduces to a generic effect
#'   with a precision matrix
#' @describeIn renderpf90 Compute the parameters of a progsf90 representation of
#'   an AR effect.
#' @export
renderpf90.ar <- function(x) {
  
  ans <- renderpf90.generic(x)
  ans$file_name = 'ar'
  
  return(ans)
}


#' Represent a symmetric matrix in triplet format
#' 
#' It only gives the lower triangular elements, and **do not** check for
#' symmetry.
#' 
#' @param x matrix.
#' @importFrom methods as
as.triplet <- function(x) {
  xsp <- Matrix::tril(as(x, 'TsparseMatrix'))
  # Note: The Matrix package counts rows and columns starting from zero
  # Thus, I add 1 to the corresponding columns
  cbind(xsp@i + 1, xsp@j + 1, xsp@x)
}