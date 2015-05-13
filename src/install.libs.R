
source_base <- "https://raw.githubusercontent.com/famuvie/breedR/master/inst/bin"
execs <- c("remlf90", "airemlf90")

if (WINDOWS) {
  source_platform <- "windows"
  execs <- paste0(execs, ".exe")
}

if (.Platform$OS.type == "unix") {
  if (is.na((file.info("/Library")$isdir && file.info("/Applications")$isdir))) {
    source_platform <- "linux"
  } else {
    source_platform <- "mac"
  }
}

f.url <- file.path(source_base, source_platform, execs)

## Download files, set execution permissions and move into place
td <- tempdir()
wd <- setwd(td)
on.exit(setwd(wd))
download.file(f.url, destfile = execs, method = 'libcurl', mode = 'wb')
Sys.chmod(execs, mode = '0774')

if ( any(file.exists(execs)) ) {
  dest <- file.path(R_PACKAGE_DIR,  paste0('bin', R_ARCH))
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(execs, dest, overwrite = TRUE)
}

