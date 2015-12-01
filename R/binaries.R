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
  
  bin.list <- progsf90_files(platform,
                             compressed = FALSE)
  
  check <- FALSE
  if (file.exists(path)) {
    check <- all(bin.list %in% list.files(path))
  }
  
  if (!check && !quiet) {
    message('Binary dependencies missing.',
            '\nWould you like to install them?\t')
    if (menu(c("Yes", "No")) == 1) {
      install_progsf90(dest = path, platform = platform)
      check <- check_progsf90(path, platform, quiet)
    }
  }
  
  return(invisible(check))
}

#' Install PROGSF90 binary dependencies
#' 
#' Copy the  binaries for the specified platform into a directory.
#' 
#' The url can be either of form http:// or of form file:// for local urls.
#' @param url where to download the files from
#' @param dest destination directory for the binaries. Default is 'bin' under
#'   the current installation dir.
#' @param platform what version of the binaries are to be installed. Default is
#'   current.
#' @param arch 
#' @param quiet logical. Whether not to display messages.
#' @export
install_progsf90 <- function(
  url      = breedr_progsf90_repo(),
  dest     = system.file('bin', package = 'breedR'),
  platform = breedR.os.type(),
  arch     = breedR.os.32or64bit(),
  quiet    = !interactive()
) {
  
  ## Check connection if URL is http:
  if (grepl("^http\\:", url) && !breedR_online()) return(FALSE)
  
  ## Binary files for this platform (packed and compressed)
  execs <- progsf90_files(platform,
                          compressed = TRUE)
  
  ## full URL for this platform and architecture
  f.url <- file.path(url, platform)
  if (platform == 'linux') # further specify arch
    f.url <- file.path(f.url, paste0(arch, 'bit'))
  
  ## Retrieve each exec to dest
  res <- sapply(execs, 
                retrieve_bin, 
                url = f.url,
                dest = dest)
  
  return(res)
}


## Download files creating dest dir if necessary
## uncompresses if necessary
## and set execution permissions
retrieve_bin <- function(f, url, dest) {
  destf <- file.path(dest, f)
  # dir.exists() does not exist in windows
  # file.exists checks for dirs as well
  if (!file.exists(dest))
    dir.create(dest, recursive = TRUE)
  
  if(grepl("^file://", url)) {
    url <- gsub("^file://", "", url)
    out <- tryCatch(
      file.copy(file.path(url, f), destf, overwrite = TRUE)
    )
    
  } else {
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
  }
  
  ## Connection issues
  if (inherits(out, 'error')) {
    unlink(destf)  # remove residual 0-byte file
    return(FALSE)
  }
  
  if (grepl("\\.zip$", f)) {
    ok <- all(unzip(destf, exdir = dest) > 0)
  } else {
    finalfs <- untar(destf, list = TRUE)
    ok <- !untar(destf, exdir = dest)
    if (ok) Sys.chmod(file.path(dest, finalfs), mode = '0744')
  }
  unlink(destf)
  
  return(ok)
}


## Return the file names of the breedR dependencies on PROGSF90 binaries 
## according to the platform
progsf90_files <- function(os = breedR.os.type(),
                           compressed = FALSE) {
  
  if (compressed) {
    ans <- paste('pf90',
                 switch(os,
                        windows = 'zip',
                        'tar.gz'),
                 sep = '.')
  } else {
    ans <- c("remlf90", "airemlf90")
    if (os == 'windows') {
      ## Ship also required dll
      ans <- c(paste0(ans, ".exe"),
               "libiomp5md.dll")
    }
  }
  
  return(ans)
}


## Check whether there is internet connection
breedR_online <- function() {
  tf <- tempfile()
  !inherits(
    suppressWarnings(
      try(download.file('http://famuvie.github.io/breedR/', tf, quiet = TRUE))
    ), 
    'try-error'
  )
}

#' Default repository for PROGSF90 binaries
breedr_progsf90_repo <- function() {
  if (!nchar(url <- Sys.getenv("PROGSF90_URL"))) {
    url <- "http://famuvie.github.io/breedR/bin"
  }
  return(url)
}
