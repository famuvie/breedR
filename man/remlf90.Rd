\name{remlf90}
\alias{remlf90}
\title{Inference with REMLF90}
\usage{
remlf90(fixed, random = NULL, var.ini = NULL, genetic = NULL,
  spatial = NULL, data, method = c("ai", "em"), debug = FALSE)
}
\arguments{
  \item{fixed}{an object of class \link{formula} (or one
  that can be coerced to that class): a symbolic
  description of the fixed effects of the model to be
  fitted. The details of model specification are given
  under 'Details'.}

  \item{random}{if not \code{NULL}, an object of class
  \link{formula} with the unstructured random effects.}

  \item{genetic}{if not \code{NULL}, a list with relevant
  parameters for an additive genetic effect; see
  'Details'.}

  \item{spatial}{if not \code{NULL}, a list with relevant
  parameters for a spatial random effect; see 'Details'.}

  \item{data}{a data frame with variables and observations}

  \item{method}{either 'ai' or 'em' for Average-Information
  or Expectation-Maximization REML respectively}
}
\value{
An object of class 'remlf90' that can be further questioned
by \code{\link{fixef}}, \code{\link{ranef}},
\code{\link{fitted}}, etc.
}
\description{
Fits a Linear Mixed Model by Restricted Maximum Likelihood
}
\details{
If either \code{genetic} or \code{spatial} are not
\code{NULL}, the model residuals are assumed to have an
additive genetic effects or a spatially structured random
effect, respectively. In those cases, code{genetic} and
\code{spatial} must be lists with named relevant
parameters.

Initial variance components can also be specified through
an additional argument \code{var.ini}. You can either use
default initial values for the variance components (see
\code{?breedR.options}) or specify custom values for
\emph{each} and \emph{all} variance components in the
model. In this case, \code{var.ini} must be a named list
with one element for each term in \code{random} with
matching names, plus one last element named \code{residual}
for the initial residual variance. Furthermore if there are
\code{genetic} or \code{spatial} effects, they must as well
include a numeric element \code{var.ini} with the initial
variance component specification for the corresponding
effect.

The available models for the genetic effect are
\code{add_animal}. \code{add_animal} stands for an additive
animal model with a given pedigree.

The available models for the spatial effect are
\code{Cappa07} and \code{AR1}. \code{Cappa07} uses a
two-dimensional tensor product of B-splines to represent
the smooth spatially structured effect. \code{AR1} uses a
kronecker product of autoregressive models for the rows and
columns.

An intercept is automatically introduced in the model
provided the user doesn't explicitly prevents it by using
\code{0} or \code{-1} in the \code{fixed} formula (as
conventional in \code{R}), \emph{and} there are no other
categorical covariates in \code{fixed}. The latter
condition is actually a limitation of (ai)remlf90 backends,
which would in any case return an estimate for each level
of the categorical covariates while returning 0 for the
intercept. It does not allow alternative parameterizations.

AI-REML is usually faster than EM-REML, and it provides
more results. Namely, standard errors of the variance
components estimates, and covariances as well. On the other
hand, is less robust than EM-REML and it usually gives
extreme results when used with the splines spatial model
(as in \code{spatial = list(model ='Cappa07')}).
}
\references{
progsf90 wiki page:
\url{http://nce.ads.uga.edu/wiki/doku.php}

E. P. Cappa and R. J. C. Cantet (2007). Bayesian estimation
of a surface to account for a spatial trend using penalized
splines in an individual-tree mixed model. \emph{Canadian
Journal of Forest Research} \strong{37}(12):2677-2688.
}
\seealso{
\code{\link[pedigreemm]{pedigree}}
}

