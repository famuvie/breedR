\name{determine.n.knots}
\alias{determine.n.knots}
\title{Determine a sensible number of knots}
\usage{
  determine.n.knots(n, cutoff = 4, rate = 0.2)
}
\arguments{
  \item{n}{An integer vector of sample sizes}
}
\value{
  An integer vector with the number of knots for each
  sample size
}
\description{
  This function computes a reasonable number of inner knots
  to be used for a basis of unidimensional B-splines
}
\references{
  Ruppert, D. (2002). Selecting the number of knots for
  penalized splines. \emph{Journal of Computational and
  Graphical Statistics} 11, 735–757.
}

