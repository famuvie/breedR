\name{as.data.frame.metagene}
\alias{as.data.frame.metagene}
\title{Coerce to a data.frame}
\usage{
  \method{as.data.frame}{metagene} (x, ...,
    exclude.founders = TRUE)
}
\arguments{
  \item{x}{a metagene object}

  \item{...}{not used}

  \item{exclude.founders}{logical: should the data.frame
  contain the genetic values of the founders?}
}
\value{
  returns a data frame with one row per individual, with
  the spatial coordinates if applicable, the pedigree
  information, the generation, the true breeding value, the
  phenotype, the sex, the spatially structured component of
  the phenotype and other internal metagene variables.
}
\description{
  This function returns the data frame with the pedigree
  and phenotypes of the simulated individuals.
}

