### Functions to handle binary dependencies  ###

#' Checks installation of binaries
#' 
#' Checks whether the binary dependencies are installed in the right directory.
#' If not, allows calling the installer
#' 
#' @param path directory to check for the presence of binaries. Default is defined in the package options, and it depends on the platform.
#' @param silent if TRUE, it won't ask whether to install missing binaries.
#' @export
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
#' 
#' @param path destination directory for the binaries. Default is taken from
#'   package option.
#' @param os.type what version of the binaries are to be installed.Default is
#'   current.
#' @export
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
