#' @importFrom stats model.matrix
#' @importMethodsFrom Matrix coerce

#' @export 
model.matrix.effect_group(x) <- function(object) {
  
  ## get the incidence matrices for all the subeffects
  mml <- lapply(object$effects, model.matrix.breedr_effect)
  
  ## confirm they have all the same number of rows
  stopifnot(length(unique(vapply(mml, nrow, 0))) == 1)
  
  ## bind by columns
  mm <- do.call(cbind, mml)
  
  return(mm)
}

#' @export 
model.matrix.breedr_effect <- function(object, ...) {
  object$incidence.matrix
}


#' @export
model.matrix.remlf90 <- function (object, ...) {
  
  ## mm and mf for fixed and diagonal effects only
  mf <- object$mf
  mm.fd <- object$mm
  
  # terms in the formula that are fixed
  fixterm.bol <- attr(attr(mf, 'terms'), 'term.types') == 'fixed'
  
  # columns in the mm corresponding to fixed and diagonal terms
  fixcol.bol <- attr(mm.fd, 'assign') %in% which(fixterm.bol)
  diacol.bol <- attr(mm.fd, 'assign') %in% which(!fixterm.bol)
  stopifnot(all(xor(fixcol.bol, diacol.bol)))  # check they are complementary
  
  
  # preallocate
  fixed  <- random <- NULL
  
  
  ## mm for fixed effects and corresponding attributes
  if( any(fixcol.bol) ) {
    fixed <- mm.fd[, fixcol.bol, drop = FALSE]
    attr(fixed, 'assign') <- attr(mm.fd, 'assign')[fixcol.bol]
    ff.bol <- names(attr(mm.fd, 'contrasts')) %in% names(which(fixterm.bol))
    if( any(ff.bol) )
      attr(fixed, 'contrasts') <- attr(mm.fd, 'contrasts')[ff.bol]
  }
  
  ## mm for diagonal effects and corresponding attributes
  for ( de.idx in which(!fixterm.bol) ) {
    de.nm <- names(fixterm.bol)[de.idx]
    decol.bol <- attr(mm.fd, 'assign') == de.idx
    stopifnot( any(decol.bol) )
    
    random[[de.nm]] <- structure(mm.fd[, decol.bol, drop = FALSE],
                                 assign = attr(mm.fd, 'assign')[decol.bol])
    
    ## Here, df.bol is TRUE at most in one element
    df.bol <- names(attr(mm.fd, 'contrasts')) %in% de.nm 
    if( any(df.bol) ) {
      cnt <- attr(mm.fd, 'contrasts')[[which(df.bol)]]
      attr(random[[de.nm]], 'contrasts') <- cnt
    }
  }
  
  
  
  ## mm for genetic effects
  if( object$components$pedigree ) {
    # Indices (in ranef) of genetic-related effects (direct and/or competition)
    gen.idx <- grep('genetic', names(object$ranef))
    
    # Incidence vector for the direct effect
    # First column of incidence matrix (and only, if model = add_animal)
    Z.direct <- as(object$effects$genetic$gen$B[,1],
                   'indMatrix')
    
    random <- c(random,
                structure(list(Z.direct),
                          names = names(object$ranef)[gen.idx[1]]))
    
    if( length(gen.idx) > 1 ) {
      
      ## model = 'competition'
      stopifnot(length(gen.idx) == 2)
      
      # Incidence matrix of competition effect is in short 16-column format
      # needs to be converted
      Z.comp <- matrix.short16(object$effects$genetic$gen$B[, 1+1:16])
      
      random <- c(random,
                  structure(list(Z.comp),
                            names = names(object$ranef)[gen.idx[2]]))
      
      ## Optional pec effect
      ## shares the same incidende matrix as the genetic competition effect
      if( exists('pec', object$ranef) ) {
        random$pec <- Z.comp
      }
      
    }
    
  }
  
  ## mm for spatial effects
  if( object$components$spatial ) {
    
    Z <- object$effects$spatial$sp$B
    
    if( !is.matrix(Z) ) {  # case AR or blocks
      ## The number of columns of the incidence matrix must be 
      ## taken from the size of the random effect
      nc <- max(object$effects$spatial$sp$U[,1])
      if( max(Z) > nc)
        stop('Incompatible dimensions between the incidence and covariance matrices in the spatial effect.')
      
      Z <- as(list(Z, nc), 'indMatrix')
    }
    random$spatial <- Z
  }
  
  
  return(list(fixed  = fixed,
              random = random))
}

