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
    warning("download failed")
    unlink(destf)  # remove residual 0-byte file
    return(FALSE)
  }

  ## unzip with internal method gave problems with R CMD check
  ## under linux. However, I want to use internal in windows.
  ## Hence, borrow solution from devtools.
  path <- decompress(destf, dest)
  unlink(destf)
  
  # Whatch out! coercion as in path > 0
  # fails when R CMD check and path starts with a slash
  return(nchar(path)>0)
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


## Borrowed from devtools
decompress <- function(src, target) {
  stopifnot(file.exists(src))
  
  if (grepl("\\.zip$", src)) {
    my_unzip(src, target)
    outdir <- getrootdir(as.vector(utils::unzip(src, list = TRUE)$Name))
    
  } else if (grepl("\\.tar$", src)) {
    utils::untar(src, exdir = target)
    outdir <- getrootdir(utils::untar(src, list = TRUE))
    
  } else if (grepl("\\.(tar\\.gz|tgz)$", src)) {
    utils::untar(src, exdir = target, compressed = "gzip")
    outdir <- getrootdir(utils::untar(src, compressed = "gzip", list = TRUE))
    
  } else if (grepl("\\.(tar\\.bz2|tbz)$", src)) {
    utils::untar(src, exdir = target, compressed = "bzip2")
    outdir <- getrootdir(utils::untar(src, compressed = "bzip2", list = TRUE))
    
  } else {
    ext <- gsub("^[^.]*\\.", "", src)
    stop("Don't know how to decompress files with extension ", ext,
         call. = FALSE)
  }
  
  file.path(target, outdir)
}


# Given a list of files, returns the root (the topmost folder)
# getrootdir(c("path/to/file", "path/to/other/thing")) returns "path/to"
getrootdir <- function(file_list) {
  slashes <- nchar(gsub("[^/]", "", file_list))
  if (min(slashes) == 0) return("")
  
  getdir(file_list[which.min(slashes)])
}

# Returns everything before the last slash in a filename
# getdir("path/to/file") returns "path/to"
# getdir("path/to/dir/") returns "path/to/dir"
getdir <- function(path)  sub("/[^/]*$", "", path)


## Adapted from devtools
my_unzip <- function(src, target, unzip = getOption("unzip")) {
  if (unzip == "internal") {
    return(utils::unzip(src, exdir = target))
  }
  
  args <- paste(
    "-oq", shQuote(src),
    "-d", shQuote(target)
  )

  ## The following is a stripped version of devtools::system_check
  ## I can't use other functions from breedR, because install.libs.R
  ## sources only os.R and binaries.R
  full <- paste(shQuote(unzip), " ", paste(args, collapse = " "), sep = "")
  result <- suppressWarnings(system(full, intern = TRUE, ignore.stderr = TRUE))
  if(is.null(status <- attr(result, "status"))) status <- 0
  if (!identical(as.character(status), "0")) {
    stop("Command failed (", status, ")", call. = FALSE)
  }
  invisible(TRUE)
}
