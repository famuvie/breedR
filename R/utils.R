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

#' Checks installation of binaries
#' 
#' Checks whether the binary dependencies are installed in the right directory.
#' If not, allows calling the installer
breedR.check.bin <- function(path = breedR.getOption('breedR.bin'),
                             silent = !interactive() ) {
  bin.list <- c('airemlf90',
                'remlf90')
  
  dir.exists <- file.exists(path)
  check <- FALSE
  if(dir.exists) {
    current.files <- sapply(strsplit(list.files(path), split='\\.'),
                            function(x) x[1])
    check <- all(!is.na(match(bin.list, current.files)))
  }
  
  if( !check && !silent ) {
    ans <- readline('Binary dependencies missing.\nWould you like to install them?\t')
    yes <- tolower(substr(ans, 1, 1) == 'y')
    
    if( yes ) {
      breedR.install.bin(path)
      check <- breedR.check.bin(path)
    }
  }
  
  return(check)
}

#' Installs binary dependencies
#' 
#' Copy the  binaries for the specified platform into a directory.
#' Default platform is current. Default path is taken from package option.
breedR.install.bin <- function(path = breedR.getOption('breedR.bin'),
                               os.type = breedR.os.type()) {
  # The repository where we have fixed and tested versions of binaries
  breedR.bin.repo <- file.path('~', 't4f', 'src', 'breedR', 'inst', 'bin')

  src.path <- file.path(breedR.bin.repo, os.type)
  bin.files <- list.files(src.path, full.names = TRUE)
  
  # Create path if necessary
  if( !file.exists(path) ) dir.create(path, recursive = TRUE)
  
  # Copy binaries into path
  out <- file.copy(from = bin.files,
                   to   = path,
                   recursive = FALSE)
  if( any(!out) ) warning('Failed to copy some of the binaries\n')
  
  return(invisible(out))
} 