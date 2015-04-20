# sp::coordinates() is an S4 function
# Register the S3 classes from breedR as S4 classes
setOldClass('breedR')
setOldClass('effect_group')
setOldClass(c('splines', 'spatial'))

#' @importFrom methods setOldClass setMethod
#' @export
setMethod('coordinates', signature = 'breedR', 
          function(obj, ...) {
            err_msg <- "This breedR object has no spatial structure.\n"
            # if( !obj$components$spatial ) {
            ## Watch out. It can be not spatial, but have coordinates
            ## assigned with coordinates <- 
            if( is.null(obj$effects$spatial$sp$coord) ) {
              ## If there is no explicit spatial structure
              ## we can check for coordinates in the genetic effect
              ## in case it was a competition model, which requires coords.
              if( !obj$components$pedigree ) {
                message(err_msg)
                return(invisible(NULL))
              } else {
                if( exists('coord', obj$effects$genetic$gen) ) {
                  return(obj$effects$genetic$gen$coord)
                } else {
                  message(err_msg)
                  return(invisible(NULL))
                }
              }
            }
            return(obj$effects$spatial$sp$coord)
          }
)

#' @export
setMethod('coordinates<-', signature = 'breedR', 
          function(object, value) {
            # sp::coordinates() performs some sanity checks
            # like coords to be positive integers, not NA, etc.
            cc <- as.data.frame(sp::coordinates(value))
            if( !object$components$spatial ) {
              object$effects$spatial <- list(name = 'none',
                                             sp = list(coord = cc))
              # Now it is a "spatial" object, although there is no
              # spatial model.
              # Not a good idea. It then confuses other methods that assume
              # that if it is 'spatial' then it has further things like
              # an incidence matrix (e.g. model.matrix)
              # object$components$spatial = TRUE
            } else {
              object$effects$spatial$sp$coord = cc
            }
            return(object)
          }
)


setMethod('coordinates', signature = 'effect_group',
          function(obj, ...) {
            
            ## Spatial effects within this group
            sp.idx <- vapply(obj$effects, inherits, TRUE, 'spatial')
            
            ## Coordinates of each one
            coord.lst <- lapply(obj$effects[sp.idx], coordinates)
            
            if (length(coord.lst) == 1)
              coord.lst <- coord.lst[[1]]
            
            return(coord.lst)
          }
)

setMethod('coordinates', signature = 'spatial',
          function(obj, ...) {
            obj$coordinates
          }
)
