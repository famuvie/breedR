\name{as.data.frame.pedigree}
\alias{as.data.frame.pedigree}
\title{Coerce to a data.frame}
\usage{
\method{as.data.frame}{pedigree}(x, ...)
}
\arguments{
  \item{x}{a \code{\link[pedigree]{pedigree}} object}

  \item{...}{not used}
}
\value{
returns a data frame with one row per individual, the first
column being the identification code, and the other two
columns are dad and mum codes respectively.
}
\description{
This function returns a \code{\link[pedigree]{pedigree}}
object in a data frame
}

