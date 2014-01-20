\name{remlf90}
\alias{remlf90}
\title{Inference with REMLF90}
\usage{
  remlf90(fixed, genetic = NULL, spatial = NULL,
    random = NULL, data, method = c("ai", "em"))
}
\arguments{
  \item{fixed}{an object of class \link{formula} (or one
  that can be coerced to that class): a symbolic
  description of the fixed effects of the model to be
  fitted. The details of model specification are given
  under 'Details'.}

  \item{genetic}{if not \code{NULL}, a list with relevant
  parameters for an additive genetic effect; see
  'Details'.}

  \item{spatial}{if not \code{NULL}, a list with relevant
  parameters for a spatial random effect; see 'Details'.}

  \item{random}{if not \code{NULL}, an object of class
  \link{formula} with the unstructured random effects.}

  \item{data}{}

  \item{method}{}
}
\value{
  An object of class 'remlf90'
}
\description{
  Fits a Linear Mixed Model by Restricted Maximum
  Likelihood
}
\details{
  If either \code{genetic} and/or \code{param} are not
  \code{NULL}, the model residuals are assumed to have an
  additive genetic effects and a spatially structured
  random effect, respectively. The relevant parameters are
  \code{model} and \code{var.ini} in both cases, and
  \code{pedigree} in the case of a genetic effect.

  The available models for the genetic effect are
  \code{add_animal}. \code{add_animal} stands for an
  additive animal model with a given pedigree.

  The available models for the spatial effect are
  \code{Cappa07}. \code{Cappa07} uses a two-dimensional
  tensor product of B-splines to represent the smooth
  spatially structured effect.
}
\references{
  progsf90 wiki page:
  \url{http://nce.ads.uga.edu/wiki/doku.php}

  E. P. Cappa and R. J. C. Cantet (2007). Bayesian
  estimation of a surface to account for a spatial trend
  using penalized splines in an individual-tree mixed
  model. \emph{Canadian Journal of Forest Research}
  \strong{37}(12):2677-2688.
}
\seealso{
  \code{\link[pedigreemm]{pedigree}}
}

