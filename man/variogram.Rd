\name{variogram}
\alias{variogram}
\title{Empirical variograms of residuals}
\usage{
variogram(x, plot = c("all", "isotropic", "anisotropic", "perspective",
  "heat", "none"), R, coord, z)
}
\arguments{
  \item{x}{a \code{breedR} result.}

  \item{plot}{character. What type of variogram is to be
  plotted. The default is 'all'. Other options are
  'isotropic', 'heat', 'perspective' and 'anisotropic'. See
  Details.}

  \item{R}{numeric. Radius of the variogram}

  \item{coord}{(optional) a two-column matrix with
  coordinates of observations}

  \item{z}{(optional) a numeric vector of values to be
  represented spatially}
}
\description{
Computes isotropic and anisotropic empirical variograms
from the residuals of a breedR model.
}
\details{
An empirical variogram computes the mean squared
differences between observations separated by a vector
\code{h}, for all possible \code{h}. An \code{isotropic}
variogram assumes that the underlying process depends only
on the (euclidean) \emph{distance} between points, but not
on the orientation or direction of \code{h}. At the other
end, an \code{anisotropic} variogram assumes that the
process might depend on the orientation, but not on the
direction. Finally, \code{heat} and \code{perspective} are
different representations of a variogram which assumes that
the process depends only on the absolute distance between
rows and columns.

Unless \code{coord} or \code{z} are specified by the user,
\code{variogram} builds the variogram with the residuals of
the model fit in \code{x}. If \code{coord} or \code{z} are
specified, then the spatial coordinates or the residuals
are respectively overrided.

This function assumes that there is at most one observation
per spatial location. Otherwise, are observations measured
at different times? should a spatial-temporal variogram be
fitted?
}

