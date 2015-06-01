# sp::coordinates() is an S4 function
# Register the S3 classes from breedR as S4 classes
setOldClass('breedR')
setOldClass('effect_group')
setOldClass(c('splines', 'spatial'))
setOldClass(c('ar', 'spatial'))
setOldClass(c('blocks', 'spatial'))

#' @importFrom methods setOldClass setMethod
#' @export
setMethod('coordinates', signature = 'breedR', 
          function(obj, ...) {
            # if( !obj$components$spatial ) {
            ## Watch out. It can be not spatial, but have coordinates
            ## assigned with coordinates <- 
            
            ## effect_groups in the list of effects
            ## (all, after refactoring)
            eg.idx <- vapply(obj$effects, inherits, TRUE, 'effect_group')
            
            ## Coordinates of each one (if any)
            ## But only those with some 
            coord.lst <- lapply(obj$effects[eg.idx], coordinates)
            coord.lst <- coord.lst[vapply(coord.lst, length, 0) > 1]
            
            ## Non-refactored components:
            cl2 <- list()
            if (!is.null(obj$effects$spatial$sp$coord)) {
              cl2 <- c(cl2, list(obj$effects$spatial$sp$coord))
            }
            if (!is.null(obj$effects$genetic$gen$coord)) {
              cl2 <- c(cl2, list(obj$effects$genetic$gen$coord))
            }
            
            ## Collected coordinates
            coord.lst <- c(coord.lst, cl2)
            
            ## If none was found, give a message and return NULL
            if (length(coord.lst) == 0) {
              message("This breedR object has no spatial structure.\n")
              return(invisible(NULL))
            }
            
            ## If there are several sets of coordinates,
            ## check that they are the same
            if (length(coord.lst) > 1) {
              if (!all.equal(coord.lst[1], coord.lst[-1]))
                stop(paste('There are different coordinate sets in the model components.',
                           'I do not know which one to use.'))
            }
            
            ## Take the first element as a representative
            coord.lst <- coord.lst[[1]]
            
            return(coord.lst)
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
