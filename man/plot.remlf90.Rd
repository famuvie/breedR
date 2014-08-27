\name{plot.remlf90}
\alias{plot.remlf90}
\title{Spatial plot of a model's fit components}
\usage{
plot.remlf90(x, type = c("phenotype", "fitted", "spatial", "fullspatial",
  "residuals"), z = NULL)
}
\arguments{
  \item{x}{A \code{breedR} object}

  \item{type}{Character. Which component is to be
  represented in the map.}

  \item{z}{Optional. A numeric vector to be plotted with
  respect to the spatial coordinates. Overrides
  \code{type}.}
}
\description{
Plots the predicted values of the component effects of the
phenotype.
}

