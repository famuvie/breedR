## Functions for checking model components
## Internal - not exported

## TODO: 
##  - check_spatial
##  - check_generic
##  - make tests for them

check_genetic <- function(model = c('add_animal', 'competition'),
                          pedigree,
                          id,
                          coordinates,
                          competition_decay = 1,
                          pec = list(present = FALSE),
                          autofill = TRUE,
                          var.ini) {
  mc <- match.call()
  
  for (arg in c('model', 'pedigree', 'id', 'var.ini')) {
    if (eval(call('missing', as.name(arg))))
      stop(paste('Argument', arg, 'required in the genetic component.'))
  }
  
  ## TODO: Verify inherits(pedigree, 'pedigree')
  
  mc$model <- match.arg(model)
  
  if (mc$model == 'competition') {
    for (arg in c('coordinates')) {
      if (eval(call('missing', as.name(arg))))
        stop(paste('Argument', arg, 'required in the genetic component.'))
    }
    
    mc$coordinates <- try(as.data.frame(coordinates))
    if (inherits(mc$coordinates, 'try-error'))
      stop('Argument coordinates in the genetic component not coercible to a data.frame')
    if (!all(vapply(coordinates, is.numeric, TRUE)))
      stop('Argument coordinates in the genetic component not numeric')
    if (ncol(coordinates) != 2)
      stop('Only two dimensions admitted for coordinates in the genetic component.')
  }
  
  stopifnot(is.numeric(competition_decay))
  stopifnot(competition_decay > 0)
  
  return(as.list(mc[-1]))
}


