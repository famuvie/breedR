\name{breedR.check.bin}
\alias{breedR.check.bin}
\title{Checks installation of binaries}
\usage{
  breedR.check.bin(path = breedR.getOption("breedR.bin"),
    silent = !interactive())
}
\arguments{
  \item{path}{directory to check for the presence of
  binaries. Default is defined in the package options, and
  it depends on the platform.}

  \item{silent}{if TRUE, it won't ask whether to install
  missing binaries.}
}
\description{
  Checks whether the binary dependencies are installed in
  the right directory. If not, allows calling the installer
}

