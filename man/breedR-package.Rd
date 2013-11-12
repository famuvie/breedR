\name{breedR-package}
\alias{breedR-package}
\alias{breedR}
\docType{package}
\title{
Frequentist and Bayesian methods for breeders and quantitative genetists
}
\description{
This package provides statistical tools to build predictive models for the breeders community.
It aims to assess the genetic value of individuals under a number of situations, including spatial autocorrelation, genetic/environment interaction and competition.
It is under active development as part of the Trees4Future project, particularly developed having forest genetic trials in mind.
But can be used for animals or other situations as well.
}
\details{
\tabular{ll}{
Package: \tab breedR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-11-12\cr
License: \tab GPL-3\cr
}
The package functionality builds up on a wrapping up of Ignacy Misztal's progsf90 suite of Fortran programs.
Particularly, the function \code{reml} performs classical Restricted-Maximum Likelihood inference by interfacing Misztal's programs with several high-level options such as spatial components, etc.
The Fortran back-end allows for fast inference on rather large datasets (hundreds of thousands of individuals) with complex pedigrees.
}
\author{
Facundo Muñoz
Maintainer: Facundo Muñoz <facundo.munoz@orleans.inra.fr>
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
~~ simple examples of the most important functions ~~
}
