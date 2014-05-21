\name{breedR.install.bin}
\alias{breedR.install.bin}
\title{Installs binary dependencies}
\usage{
  breedR.install.bin(path = breedR.getOption("breedR.bin"),
    os.type = breedR.os.type())
}
\arguments{
  \item{path}{destination directory for the binaries.
  Default is taken from package option.}

  \item{os.type}{what version of the binaries are to be
  installed.Default is current.}
}
\description{
  Copy the binaries for the specified platform into a
  directory.
}

