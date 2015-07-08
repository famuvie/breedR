#' @importFrom stats model.matrix
#' @importMethodsFrom Matrix coerce

#' @export 
model.matrix.effect_group <- function(object, ...) {
  
  ## get the incidence matrices for all the subeffects
  mml <- lapply(object$effects, model.matrix, ...)
  
  ## confirm they have all the same number of rows
  stopifnot(length(unique(vapply(mml, nrow, 0))) == 1)
  
  ## bind by columns
  # No: keep components separated (like direct and competition)
  # mm <- do.call(Matrix::cBind, mml)
  
  return(mml)
}

#' @export 
model.matrix.splines <- function(object, fullgrid = FALSE) {

  if (!fullgrid) return(model.matrix.breedr_effect(object))
  
  coord <- coordinates(object)
  obs.loc <- loc_grid(coord, autofill = TRUE)
  grid <- expand.grid(obs.loc, KEEP.OUT.ATTRS = FALSE)
  
  inc.mat <- bispline_incidence(object$knots,
                                grid,
                                object$degree + 1,
                                sparse = TRUE)
  ans <- structure(inc.mat,
                   coordinates = grid)
}

#' @export 
model.matrix.ar <- function(object, fullgrid = FALSE) {
  
  if (!fullgrid) return(model.matrix.breedr_effect(object))
  
  grid.def <- attr(object, 'grid')
  pos.idx <- lapply(lapply(grid.def$length, seq_len), `-`, 1)
  pos.m <- mapply(`*`, lapply(grid.def$step, identity), pos.idx,
                  SIMPLIFY = FALSE)
  pos.loc <- mapply(`+`, lapply(grid.def$origin, identity), pos.m,
                    SIMPLIFY = FALSE)
  grid <- expand.grid(pos.loc)
  
  inc.mat <- as(seq_len(nrow(grid)), 'indMatrix')
  
  ans <- structure(inc.mat,
                   coordinates = grid)
}


#' @export 
model.matrix.blocks <- function(object, fullgrid = FALSE) {

  ## I don't really have a way to infer the corresponding block
  ## for unobserved locations.
  ans <- model.matrix.breedr_effect(object)
  
  if (fullgrid)
    ans <- structure(ans,
                     coordinates = coordinates(object))
  
  return(ans)
}

#' @export 
model.matrix.splines <- function(object, fullgrid = FALSE) {

  if (!fullgrid) return(model.matrix.breedr_effect(object))
  
  coord <- coordinates(object)
  obs.loc <- loc_grid(coord, autofill = TRUE)
  grid <- expand.grid(obs.loc, KEEP.OUT.ATTRS = FALSE)
  
  inc.mat <- bispline_incidence(object$knots,
                                grid,
                                object$degree + 1,
                                sparse = TRUE)
  ans <- structure(inc.mat,
                   coordinates = grid)
}


#' @export 
model.matrix.breedr_effect <- function(object, ...) {

  ans <-  object$incidence.matrix
  
  ## Ensure a matrix as a result
  if (!inherits(ans, 'matrix') && !inherits(ans, 'Matrix'))
    ans <- model.matrix(~x-1, data.frame(x = ans))
  
  return(ans)
}


#' @export
#' @family extraction functions
model.matrix.remlf90 <- function (object, ...) {
  
  mml <- lapply(object$effects, model.matrix)
  
  ## Fixed effects are matrices, while effect_groups
  ## are lists of matrices (even of length 1)
  
  ## Turn fixed effects to a named list of one matrix
  fix.idx <- sapply(object$effects, effect_type) == 'fixed'
  mml[fix.idx] <- lapply(which(fix.idx), function(x) mml[x])

  ## Flatten the list of lists keeping proper names  
  stopifnot(all(vapply(mml, is.list, TRUE)))
  mml <- unlist(mml, recursive = FALSE)
  names(mml) <- get_efnames(object$effects)
  
  ## confirm they have all the same number of rows
  stopifnot(length(unique(vapply(mml, nrow, 0))) == 1)
  
  
  #   ## mm and mf for fixed and diagonal effects only
  #   mf <- model.frame.remlf90(object)
  #   mm.fd <- object$mm
  #   
  #   # terms in the formula that are fixed
  #   fixterm.bol <- attr(attr(mf, 'terms'), 'term.types') == 'fixed'
  #   
  #   # columns in the mm corresponding to fixed and diagonal terms
  #   fixcol.bol <- attr(mm.fd, 'assign') %in% which(fixterm.bol)
  #   diacol.bol <- attr(mm.fd, 'assign') %in% which(!fixterm.bol)
  #   stopifnot(all(xor(fixcol.bol, diacol.bol)))  # check they are complementary
  #   
  #   
  #   # preallocate
  #   fixed  <- random <- NULL
  #   
  #   
  #   ## mm for fixed effects and corresponding attributes
  #   if( any(fixcol.bol) ) {
  #     fixed <- mm.fd[, fixcol.bol, drop = FALSE]
  #     attr(fixed, 'assign') <- attr(mm.fd, 'assign')[fixcol.bol]
  #     ff.bol <- names(attr(mm.fd, 'contrasts')) %in% names(which(fixterm.bol))
  #     if( any(ff.bol) )
  #       attr(fixed, 'contrasts') <- attr(mm.fd, 'contrasts')[ff.bol]
  #   }
  #   
  #   ## mm for diagonal effects and corresponding attributes
  #   for ( de.idx in which(!fixterm.bol) ) {
  #     de.nm <- names(fixterm.bol)[de.idx]
  #     decol.bol <- attr(mm.fd, 'assign') == de.idx
  #     stopifnot( any(decol.bol) )
  #     
  #     random[[de.nm]] <- structure(mm.fd[, decol.bol, drop = FALSE],
  #                                  assign = attr(mm.fd, 'assign')[decol.bol])
  #     
  #     ## Here, df.bol is TRUE at most in one element
  #     df.bol <- names(attr(mm.fd, 'contrasts')) %in% de.nm 
  #     if( any(df.bol) ) {
  #       cnt <- attr(mm.fd, 'contrasts')[[which(df.bol)]]
  #       attr(random[[de.nm]], 'contrasts') <- cnt
  #     }
  #   }
  #   
  #   
  #   
  #   ## mm for genetic effects
  #   if( object$components$pedigree ) {
  #     # Indices (in ranef) of genetic-related effects (direct and/or competition)
  #     gen.idx <- grep('genetic', names(object$ranef))
  #     
  #     # Incidence vector for the direct effect
  #     # First column of incidence matrix (and only, if model = add_animal)
  #     Z.direct <- as(object$effects$genetic$gen$B[,1],
  #                    'indMatrix')
  #     
  #     random <- c(random,
  #                 structure(list(Z.direct),
  #                           names = names(object$ranef)[gen.idx[1]]))
  #     
  #     if( length(gen.idx) > 1 ) {
  #       
  #       ## model = 'competition'
  #       stopifnot(length(gen.idx) == 2)
  #       
  #       # Incidence matrix of competition effect is in short 16-column format
  #       # needs to be converted
  #       Z.comp <- matrix.short16(object$effects$genetic$gen$B[, 1+1:16])
  #       
  #       random <- c(random,
  #                   structure(list(Z.comp),
  #                             names = names(object$ranef)[gen.idx[2]]))
  #       
  #       ## Optional pec effect
  #       ## shares the same incidende matrix as the genetic competition effect
  #       if( exists('pec', object$ranef) ) {
  #         random$pec <- Z.comp
  #       }
  #       
  #     }
  #     
  #   }
  #   
  #   ## mm for spatial effects
  #   ## Remove after refactoring
  #   if (object$components$spatial) {
  #     
  #     ## Only for not-yet refactored objects
  #     if (!inherits(object$effects$spatial, 'effect_group')) {
  #       
  #       Z <- object$effects$spatial$sp$B
  #       
  #       if( !is.matrix(Z) & !inherits(Z, 'Matrix')) {  # case AR or blocks
  #         ## The number of columns of the incidence matrix must be 
  #         ## taken from the size of the random effect
  #         nc <- max(object$effects$spatial$sp$U[,1])
  #         if( max(Z) > nc)
  #           stop('Incompatible dimensions between the incidence and covariance matrices in the spatial effect.')
  #         
  #         Z <- as(list(Z, nc), 'indMatrix')
  #       }
  #       random$spatial <- Z
  #     }
  #   }
  #   
  #   ## mm for refactored effects
  #   rf.idx <- vapply(object$effects, inherits, TRUE, 'effect_group')
  #   
  #   if (any(rf.idx)) {
  #     rf.mm <- lapply(object$effects[rf.idx], model.matrix.effect_group)
  #     random <- c(random,
  #                 rf.mm)
  #   }
  
  return(mml)
}

