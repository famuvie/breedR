\name{loc_grid}
\alias{loc_grid}
\title{Lattice of spatial locations}
\usage{
loc_grid(coord, autofill)
}
\arguments{
  \item{coord}{matrix, data.frame or list of row and column
  coordinates of observational units}

  \item{autofill}{locigcal, whether to perform automatic
  filling of empty rows/columns}
}
\value{
list of row and column coordinates of spatial units

This functions converts observational coordinates into
spatial coordinates.  While in a observational dataset
there can possibly be many or missing observations at one
single location, this function returns the list of row and
column coordinates of a grid which contains all the
observations
}
\description{
Returns a list row and column coordinates of observations
}

