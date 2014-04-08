\name{breedR.option}
\alias{breedR.getOption}
\alias{breedR.option}
\alias{breedR.options}
\alias{breedR.setOption}
\title{Set and get global options for breedR}
\usage{
breedR.getOption(option = c("ar.eval", "splines.nok",
  "default.initial.variance"))

breedR.setOption(...)
}
\arguments{
  \item{...}{Option and value, like \code{option=value} or
  \code{option, value}; see the Examples.}

  \item{option}{The option to get. If missing or
  \code{NULL}, then \code{breedR.getOption} will display
  the current defaults, otherwise, \code{option} must be
  one of

  \code{ar.eval}: numeric vector of values in (-1, 1) where
  the autoregressive parameters should be evaluated if not
  otherwise specified

  \code{splines.nok}: a function of the number of
  individuals in a row which gives the number of knots
  (nok) to be used for a splines model, if not otherwise
  specified

  \code{default.initial.variance}: a default value for all
  variance components}
}
\description{
Set and get global options for breedR. The options are
stored in the variable \code{breedR.options} in the
\code{.GlobalEnv}-environment.
}
\examples{
## Set default values for the autoregressive parameters
breedR.setOption("ar.eval", 3*(-3:3)/10)
## alternative format
breedR.setOption(ar.eval = 3*(-3:3)/10)
## check it
breedR.getOption("ar.eval")
}

