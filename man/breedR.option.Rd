\name{breedR.option}
\alias{Set}
\alias{The}
\alias{\code{.GlobalEnv}-environment.}
\alias{\code{breedR.options}}
\alias{and}
\alias{are}
\alias{breedR.}
\alias{breedR.getOption}
\alias{breedR.option}
\alias{breedR.options}
\alias{breedR.setOption}
\alias{for}
\alias{get}
\alias{global}
\alias{in}
\alias{options}
\alias{stored}
\alias{the}
\alias{variable}
\alias{}}
\title{Set and get global options for breedR}
\usage{
breedR.getOption(option = c("ar.eval", "splines.nok"))
}
\arguments{
  \item{...}{Option and value, like \code{option=value} or
  \code{option, value}; see the Examples.}

  \item{option}{The option to get. If \code{option = NULL}
  then \code{breedR.getOption} will display the current
  defaults, otherwise, \code{option} must be one of

  ar.eval: numeric vector of values in (-1, 1) where the
  autoregressive parameters should be evaluated if not
  otherwise specified

  splines.nok: a function of the number of individuals in a
  row which gives the number of knots (nok) to be used for
  a splines model, if not otherwise specified}
}
\description{
Set and get global options for breedR
}
\examples{
## Set default values for the autoregressive parameters
breedR.setOption("ar.eval", 3*(-3:3)/10)
## alternative format
breedR.setOption(ar.eval = 3*(-3:3)/10)
## check it
breedR.getOption("ar.eval") }
}

