### Functions to download and install binary dependencies  ###

#' Checks installation of PROGSF90 binaries
#' 
#' Checks whether the binary dependencies are installed in the right directory. 
#' If not, allows calling the installer
#' 
#' This function does not check whether the binaries are for the right platform 
#' or architecture. It only checks the presence of files with the expected 
#' names.
#' 
#' @param path directory to check for the presence of binaries. Default is
#'   defined in the package options, and it depends on the platform.
#' @param quiet if TRUE, it won't ask whether to install missing binaries.
#' @export
check_progsf90 <- function(path = breedR.getOption('breedR.bin'),
                           platform = breedR.os.type(),
                           quiet = !interactive() ) {

  bin.list <- progsf90_files(platform)
  
  check <- FALSE
  if (file.exists(path)) {
    current.files <- sapply(strsplit(list.files(path), split='\\.'),
                            function(x) x[1])
    check <- all(!is.na(match(bin.list, current.files)))
  }
  
  if (!check && !quiet) {
    message('Binary dependencies missing.\nWould you like to install them?\t')
    ans <- readline()
    yes <- tolower(substr(ans, 1, 1) == 'y')
    
    if( yes ) {
      install_progsf90(dest = path)
      check <- check_progsf90(path, quiet)
    }
  }
  
  return(invisible(check))
}

#' Install PROGSF90 binary dependencies
#' 
#' Copy the  binaries for the specified platform into a directory.
#' 
#' @param url where to download the files from
#' @param dest destination directory for the binaries. Default is 'bin' under
#'   the current installation dir.
#' @param platform what version of the binaries are to be installed. Default is
#'   current.
#' @param arch 
#' @param quiet logical. Whether not to display messages.
#' @export
install_progsf90 <- function(
  url = "http://famuvie.github.io/breedR/bin",
  dest   = system.file('bin', package = 'breedR'),
  platform = breedR.os.type(),
  arch  = paste0(breedR.os.32or64bit(), 'bit'),
  quiet = !interactive()
) {
  
  execs <- progsf90_files(platform)
  
  f.url <- file.path(url, platform, arch)
  if (platform == 'mac')  # remove arch for mac
    f.url <- dirname(f.url)
  
  res <- sapply(execs, 
                retrieve_bin, 
                url = f.url,
                dest = dest)
  
  return(res)
}


## Download files creating dest dir if necessary
## and set execution permissions
retrieve_bin <- function(f, url, dest) {
  destf <- file.path(dest, f)
  if (!dir.exists(dest))
    dir.create(dest, recursive = TRUE)
  out <- tryCatch(
    download.file(
      url = file.path(url, f),
      destfile = destf,
      mode = 'wb',
      cacheOK = FALSE,
      quiet = TRUE
    ),
    error = identity
  )
  
  ## Connection issues
  if (inherits(out, 'error')) {
    unlink(destf)  # remove residual 0-byte file
    return(FALSE)
  }
  
  Sys.chmod(destf, mode = '0744')
  return(destf)
}


## Return the file names of the breedR dependencies on PROGSF90 binaries 
## according to the platform
progsf90_files <- function(os = breedR.os.type()) {

  ans <- c("remlf90", "airemlf90")
  if (os == 'windows') {
    ## Ship also required dll
    ans <- c(paste0(ans, ".exe"),
             "libiomp5md.dll")
  }
  
  ans
}
