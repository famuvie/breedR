## Internal utility functions
## Not exported

#' lme4-style formulas
#' 
#' Transform the separated fixed and random formulas
#' into the single formula with lme4 syntaxis
lme4_fml <- function(fix, rnd, rm_int = TRUE) {
  rnd.terms <- attr(terms(rnd), 'term.labels')
  rnd.terms.lme4 <- paste('(1|', rnd.terms, ')', sep ='')
  int <- ifelse(rm_int, '-1', '')
  rnd.upd <- paste('~ .', int, paste('+', rnd.terms.lme4, collapse = ' '))
  fml.lme4 <- update(fix, rnd.upd)
  return(fml.lme4)
}


breedR.is.element <- function(name, alist) {
  ## return TRUE if element with name NAME is a member of LIST and
  ## the value is non null and not NA.
  if (any(names(alist) == name)) {
    idx = which(names(alist) == name)
    if (!is.null(alist[[idx]]) && !suppressWarnings(is.na(alist[[idx]]))) {
      return (TRUE)
    } else {
      return (FALSE)
    }
  } else {
    return (FALSE)
  }
}

breedR.get.element <-  function(name, alist) {
  if (breedR.is.element(name, alist)) {
    return (alist[[which(names(alist) == name)]])
  } else {
    return (NULL)
  }
}


#' Fit some model
breedR.result <- function() {
  res  <- suppressWarnings(remlf90(fixed  = phe_X ~ gg,
                                   genetic = list(model = 'add_animal', 
                                                  pedigree = globulus[,1:3],
                                                  id = 'self'), 
                                   spatial = list(model = 'AR', 
                                                  coord = globulus[, c('x','y')],
                                                  rho = c(.85, .8)), 
                                   data = globulus))
  return(res)
}


#' Return platform string
#' 
#' Return whether the OS is either \code{windows}, \code{linux} or \code{mac}
#' Inspired in INLA's os.R functions
breedR.os.type <- function() {
  
  type <- .Platform$OS.type

  if (type == "windows") {
    os <- type
  } else if (type == "unix") {
    mac.dirs <- file.info("/Library")$isdir && file.info("/Applications")$isdir
    os <- ifelse(is.na(mac.dirs), 'linux', 'mac')
  } else os <- 'else'
    
  return(os)
}