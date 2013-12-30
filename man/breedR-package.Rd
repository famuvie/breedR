\docType{package}
\name{breedR-package}
\alias{breedR}
\alias{breedR-package}
\title{Frequentist and Bayesian methods for breeders, quantitative genetists and
forest genetic resources analysts.}
\description{
  This package provides statistical tools to build
  predictive models for the breeders, quantitative
  genetists and forest genetic resources analysts
  communities.  It aims to assess the genetic value of
  individuals under a number of situations, including
  spatial autocorrelation, genetic/environment interaction
  and competition. It is under active development as part
  of the Trees4Future project, particularly developed
  having forest genetic trials in mind. But can be used for
  animals or other situations as well.
}
\details{
  The package functionality builds up on a wrapping up of
  Ignacy Misztal's progsf90 suite of Fortran programs.
  Particularly, the function \code{reml} performs classical
  Restricted-Maximum Likelihood inference by interfacing
  Misztal's programs with several high-level options such
  as spatial components, etc. The Fortran back-end allows
  for fast inference on rather large datasets (hundreds of
  thousands of individuals) with complex pedigrees.
}
\examples{
# Load, summarize and visualize data
  data(m4)
  summary(m4)
  plot(m4)

  # Fit Mixed Model using REML
  res.f90 <- remlf90(formula = phe_X ~ gen,
                     genetic = list(model = 'add_animal',
                                    var.ini = 10,
                                    pedigree = pedigree(m4),
                                    id = 'self'),
                     data = as.data.frame(m4))

  # Summary of results
  summary(res.f90)

  # Observed phenotypes vs. Fitted values
  qplot(phe_X, fitted(res.f90), color=gen, data = as.data.frame(m4)) + geom_abline(int=0, slope=1)
}
\references{
  Most functionality in the package is based on Ignacy
  Misztal's suite of Fortran programs for mixed model
  computations in breeding.
  \url{http://nce.ads.uga.edu/wiki/doku.php}

  M. Lynch \& B. Walsh (1998). \emph{Genetics and Analysis
  of Quantitative Traits}. Sinauer Associates, Inc.
}
\seealso{
  \code{\link[pedigreemm]{pedigreemm}}
}
\keyword{package}

