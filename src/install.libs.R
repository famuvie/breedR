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

install_progsf90(dest = file.path(R_PACKAGE_DIR, 'bin'))
