# Script which will be run as part of the installation used to automate install 
# of PROGSF90 binaries The script is run in a separate R environment containing 
# the following variables: R_PACKAGE_NAME (the name of the package), 
# R_PACKAGE_SOURCE (the path to the source directory of the package), 
# R_PACKAGE_DIR (the path of the target installation directory of the package), 
# R_ARCH (the arch-dependent part of the path, often empty), SHLIB_EXT (the
# extension of shared objects) and WINDOWS (TRUE on Windows, FALSE elsewhere). 
# REF:
# http://cran.univ-paris1.fr/doc/manuals/r-release/R-exts.html#Package-subdirectories


source('../R/binaries.R')
source('../R/os.R')

# cat("PROGSF90_URL:", Sys.getenv("PROGSF90_URL"))
message("Downloading PROGSF90 from:\n", breedr_progsf90_repo())
install_progsf90(dest = file.path(R_PACKAGE_DIR, 'bin'))
if (WINDOWS) {
  ## This is run on the builder machine. I don't know the user's architecture
  ## at this point. Nor do I know whether there will be Internet connection
  ## available at load time. Therefore, I pack both versions of PROGSF90.
  install_progsf90(dest = file.path(R_PACKAGE_DIR, 'bin', '32bit'), arch = 32)
  install_progsf90(dest = file.path(R_PACKAGE_DIR, 'bin', '64bit'), arch = 64)
  
} else {
  install_progsf90(dest = file.path(R_PACKAGE_DIR, 'bin'))
}
