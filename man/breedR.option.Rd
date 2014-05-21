\name{breedR.option}
\alias{breedR.getOption}
\alias{breedR.option}
\alias{breedR.options}
\alias{breedR.setOption}
\title{Set and get global options for breedR}
\usage{
  breedR.getOption(option = c("ar.eval", "splines.nok", "default.initial.variance", "col.seq", "col.div", "breedR.bin"))

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
  variance components

  \code{col.seq}: a vector with the specification of
  default extreme breedR col for sequential scales in
  spatial quantitative plots. See Details.

  \code{col.div}: a vector with the specification of
  default extreme breedR col for diverging scales in
  spatial quantitative plots. See Details.

  \code{breedR.bin}: full path for breedR backend binaries}
}
\description{
Set and get global options for breedR. The options are
stored in the variable \code{breedR.options} in the
\code{.GlobalEnv}-environment.
}
\details{
Sequential scales are used for variables not necessarily
centered such as a response variable, or the fitted values
of a model. The colour scale is built as a gradient between
two extreme colours which are specified as hex codes or
colour names in the option \code{col.seq}.

Diverging scales are used for plots such as residuals,
centered (hopefully) around zero, with positive and
negative values represented with different colours whose
intensity is linked to the magnitude. The option
\code{col.div} is a vector of two hex codes or colour names
of the most intense colours.
}
\examples{
## Set default values for the autoregressive parameters
breedR.setOption("ar.eval", 3*(-3:3)/10)
## alternative format
breedR.setOption(ar.eval = 3*(-3:3)/10)
## check it
breedR.getOption("ar.eval")
}

