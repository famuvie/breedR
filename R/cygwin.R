## Nothing to export

# Return TRUE if CYGWIN seems installed in the given dir
breedR.cygwin.check <- function(path = breedR.getOption("cygwin")) {
  return (file.exists(path) && file.info(path)$isdir)
}


# Check whether cygwin/bin is in the PATH environment variable
# and inserts it if necessary. 
breedR.cygwin.setPATH <- function(path = breedR.getOption("cygwin")) {
  cygbin <- file.path(path, 'bin')
  if( !grepl(cygbin, tolower(Sys.getenv('PATH'))) ) {
    path <- paste(cygbin, Sys.getenv('PATH'), sep = .Platform$path.sep)
    Sys.setenv(PATH = path)
  }
}
