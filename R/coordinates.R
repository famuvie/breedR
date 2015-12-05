# sp::coordinates() is an S4 function
# Register the S3 classes from breedR as S4 classes
setOldClass('breedR')
setOldClass('remlf90')
setOldClass('effect_group')
setOldClass(c('blocks', 'spatial'))
setOldClass(c('ar', 'spatial'))
setOldClass(c('splines', 'spatial'))
setOldClass(c('additive_genetic_competition',
              'additive_genetic', 'genetic',
              'competition', 'spatial'))
setOldClass(c("permanent_environmental_competition",
              "competition", "spatial"))

#' breedR coordinates methods
#' 
#' @param obj an object of the corresponding class
#' @param object an object of the corresponding class
#' @param value 2-column matrix or data frame with coordinates
#' @param ... not used.
#' @name coordinates_breedR
NULL

#' @importFrom methods setOldClass setMethod
#' @rdname coordinates_breedR
#' @export
setMethod('coordinates', signature = 'breedR', 
          function(obj, ...) {
            # if( !obj$components$spatial ) {
            ## Watch out. It can be not spatial, but have coordinates
            ## assigned with coordinates <- 
            
            ## effect_groups in the list of effects
            eg.idx <- which(lapply(obj$effects, effect_type) == 'random')
            
            ## Coordinates of each one (if any)
            ## But only those with some 
            coord.lst <- lapply(obj$effects[eg.idx], coordinates)
            coord.lst <- coord.lst[vapply(coord.lst, length, 0) > 1]
            
            if (length(coord.lst) == 0) {
              ## Check for an specific list item
              if ('coordinates' %in% names(obj)) {
                coord.lst <- list(obj$coordinates)
              } else {
                ## If none was found, give a message and return NULL
                message("This breedR object has no spatial structure.\n")
                return(invisible(NULL))
              }
            }
            
            ## If there are several sets of coordinates,
            ## check that they are the same
            if ((nc <- length(coord.lst)) > 1) {
              for (i in 2:nc){
                if (!isTRUE(all.equal(coord.lst[[1]], coord.lst[[i]])))
                  stop(paste('There are different coordinate sets in the model components.',
                             'I do not know which one to use.'))
              }
            }
            
            ## Take the first element as a representative
            coord.lst <- coord.lst[[1]]
            
            return(coord.lst)
          }
)

#' @rdname coordinates_breedR
#' @export
setMethod('coordinates<-', signature = 'breedR', 
          function(object, value) {
            # sp::coordinates() performs some sanity checks
            # like coords to be positive integers, not NA, etc.
            cc <- as.data.frame(sp::coordinates(value))
            
            if (is.null(suppressMessages(coordinates(object)))) {
              object$coordinates <- cc
            } else {
              ## This object already has coordinates
              ## give message and do nothing
              message("This breedR object already has coordinates.\n")
            }
            return(object)
          }
)


#' @rdname coordinates_breedR
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

#' @rdname coordinates_breedR
setMethod('coordinates', signature = 'spatial',
          function(obj, ...) {
            obj$coordinates
          }
)
