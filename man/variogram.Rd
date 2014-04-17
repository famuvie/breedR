\name{variogram}
\alias{variogram}
\title{Empirical variograms of residuals}
\usage{
variogram(x, R, plot = c("all", "isotropic", "anisotropic", "perspective",
  "heat", "none"))
}
\description{
Computes isotropic and anisotropic empirical variograms
from the residuals of a breedR model.
}
\details{
This assumes that there is at most one observation per
spatial location. Otherwise, are observations measured at
different times? should a spatial-temporal variogram be
fitted?
}

